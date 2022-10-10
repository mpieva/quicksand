#!/usr/bin/env nextflow

nextflow.enable.dsl = 1

red = "\033[0;31m"
white = "\033[0m"
cyan = "\033[0;36m"
yellow = "\033[0;33m"
standard_run = true

log.info """
[quicksand]: Execution started: ${workflow.start.format('dd.MM.yyyy HH:mm')} ${cyan}

  ============================================================
  ==========================  =============================  =
  ==    ==  =  ==  ===   ===  =  ===   ====   ===  ========  =
  =  =  ==  =  ======  =  ==    ===  =  ==  =  ==     ===    =
  =  =  ==  =  ==  ==  =====   =====  =======  ==  =  ==  =  =
  ==    ==  =  ==  ==  =====    =====  ====    ==  =  ==  =  =
  ====  ==  =  ==  ==  =  ==  =  ==  =  ==  =  ==  =  ==  =  =
  ====  ===    ==  ===   ===  =  ===   ====    ==  =  ===    =
  ============================================================                  
  ${white}${workflow.manifest.description} ${cyan}~ Version ${workflow.manifest.version} ${white}

 --------------------------------------------------------------
"""

//
//
// Functions used within the pipeline
//
//

def exit_with_error_msg(error, text){
    println "[quicksand]: ${red}${error}: ${text}${white}"
    exit 0
}
def get_warn_msg(text){
    return "[quicksand]: ${yellow}(WARN): ${text}${white}"
}
def get_info_msg(text){
    return "[quicksand]: ${text}"
}
def has_ending(file, extension){
    return extension.any{ file.toString().toLowerCase().endsWith(it) }
}
def exit_missing_required(flag){
    exit_with_error_msg("ArgumentError", "Missing required Argument ${flag}")
}

//
//
// Help
//
//

if (params.help){
    print file("$baseDir/assets/help.txt").text
    exit 0
}

//
//
// Channel assignment
//
//


// Main File Input Channels 
inbam      = params.bam            ? Channel.fromPath("${params.bam}",                   checkIfExists:true)  : Channel.empty()
indexfile  = params.rg             ? Channel.fromPath("${params.rg}",                    checkIfExists:true)  : Channel.empty()
splitdir   = params.split          ? Channel.fromPath("${params.split}/*",               checkIfExists:true)  : Channel.empty()
genomesdir = params.genomes        ? Channel.fromPath("${params.genomes}", type:'dir',   checkIfExists:true)  : Channel.empty()
database   = params.db             ? Channel.fromPath("${params.db}", type:'dir',        checkIfExists:true)  : Channel.empty()
bedfiledir = params.bedfiles       ? Channel.fromPath("${params.bedfiles}",type:'dir',   checkIfExists:true)  : Channel.empty()


// Additional File Channels
taxid = new File("${params.genomes}/taxid_map.tsv").exists() ? Channel.fromPath("${params.genomes}/taxid_map.tsv", type:'file') : Channel.fromPath("${baseDir}/assets/taxid_map_example.tsv", type:'file') 

if ( !new File("${params.genomes}/taxid_map.tsv").exists() && params.genomes) {
    log.info get_warn_msg("The file 'taxid_map.tsv' is missing in your genomes dir! Using fallback ${baseDir}/assets/taxid_map_example.tsv")
}    

// for only-fixed runs
report = new File('final_report.tsv').exists() && params.rerun ? Channel.fromPath('final_report.tsv', type:'file') : Channel.empty()
report
  .splitCsv(sep:'\t', header:true)
  .map{row -> [row.RG, row.Family, row]}
  .into{report;final_report}

outdir = params.rerun ? Channel.fromPath('out/*/1-extracted/*.bam', type: 'file') : Channel.empty()
outdir
  .map{it -> [it.baseName.split('-')[-1], it, it.baseName.split('_')[0]]}
  .set{outdir}

//
//
// Validation
//
//

if(params.taxlvl !in ['f','o']){
    exit_with_error_msg("ArgumentError","taxlvl must be one of [o, f] not ${params.taxlvl}")
}

if( params.rerun && params.fixed ){ standard_run = false }

if(params.split && (params.bam || params.rg) && standard_run){
    log.info get_info_msg("Use: nextflow run mpieva/quicksand {--rg FILE --bam FILE | --split DIR}")
    exit_with_error_msg("ArgumentError", "Too many arguments")
} 
if(!params.split && !(params.bam && params.rg) && standard_run){
    log.info get_info_msg("Use: nextflow run mpieva/quicksand {--rg FILE --bam FILE | --split DIR}")
    exit_with_error_msg("ArgumentError", "Too few arguments")
}
if(params.rerun && !(params.fixed)){
    log.info get_info_msg("Use --rerun together with --fixed")
    exit_with_error_msg("ArgumentError", "Too few arguments")
}
if(standard_run && !params.genomes){ exit_missing_required('--genomes') }
if(standard_run && !params.db){   exit_missing_required('--db')      }
if(standard_run && !params.bedfiles){ exit_missing_required('--bedfiles')}


//
//
// Pipeline
//
//

process splitBam {
    container (workflow.containerEngine ? "merszym/splitbam:v0.1.6" : null)
    publishDir 'split', mode: 'copy', pattern: '*.bam'
    publishDir 'split', mode: 'copy', pattern: '*.txt'
    label 'local'

    input:
    path 'input.bam' from inbam
    path 'indices.tsv' from indexfile

    output:
    path '*' into splitbam_out mode flatten

    when:
    params.rerun == false && params.bam && params.rg   

    script:
    """
    splitbam -s -c $params.compression_level -f indices.tsv --minscore 10 --maxnumber 0 input.bam > splittingstats.txt
    """
}

// If split is defined, start the pipeline here
if(params.split){
    splitbam_out = splitdir
} 

//create the meta-map
splitbam_out
   .map{[['RG': it.baseName, 'Reference':'best', 'ref_path':null, 'ExtractLVL':params.taxlvl], it]}      
   .set{splitbam_out}

//convert fastq-input
splitbam_out
    .branch{
        bam: it[1].getExtension() == "bam"
        fastq: has_ending(it[1], ["fastq","fastq.gz","fq","fq.gz"])
        split: it[1].name =~ /split.*stats/
        fail: true
    }
    .set{splitbam_out}

splitbam_out.fail
    .view{get_warn_msg("${it[1]} omitted. File has neither bam nor fastq-ending!")}

crosscont_in = splitbam_out.split?: null 

process EstimateCrossContamination{
    container (workflow.containerEngine ? "merszym/bam_deam:nextflow" : null)
    label 'local'

    input:
    tuple meta, 'splittingstats.txt' from crosscont_in

    output:
    tuple meta, "${meta.RG}.CC.txt" into crosscont_out 

    when:
    params.rerun == false   

    script:
    """
    cross_cont.py splittingstats.txt > \"${meta.RG}.CC.txt\"
    """
}

crosscont_out
    .filter{it[1].text != ''}
    .collectFile(storeDir: '.', newLine:true) {meta, file ->
        [ "cc_estimates.txt", ["# Cross contamination calculated from File: ${meta.RG}.txt"  ,file.text].join('\n') ]
    }

process fastq2Bam{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    tag "$meta.RG"
    label 'local'
    
    input:
    tuple meta, "${meta.RG}.fq" from splitbam_out.fastq
    
    output:
    tuple meta, "${meta.RG}.bam" into fastq2bam_out
    
    when:
    params.rerun == false   

    script:
    """
    samtools import -0 \"${meta.RG}.fq\" -o \"${meta.RG}.bam\"
    """
}

splitbam_out.bam.mix(fastq2bam_out).set{filterbam_in}

process filterBam {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    label 'local'
    tag "$meta.RG:Flag:${params.bamfilterflag}"

    input:
    tuple meta, "${meta.RG}.bam" from filterbam_in

    output:
    tuple meta, "${meta.RG}.filtered.bam", stdout, 'filtercount.txt' into filterbam_out

    when:
    params.rerun == false   

    script:
    """
    samtools view -c -F 128 \"${meta.RG}.bam\"
    samtools view -b -u -F ${params.bamfilterflag} -o \"${meta.RG}.filtered.bam\" \"${meta.RG}.bam\"
    samtools view -c \"${meta.RG}.filtered.bam\" > 'filtercount.txt'
    """
}

//add splitcounts to meta
filterbam_out
    .map{meta,bam,raw,filter -> [meta+['ReadsRaw':raw.trim(),'ReadsFiltered':filter.text.trim()], bam]}
    .set{filterlength_in}

process filterLength {
    container (workflow.containerEngine ? "merszym/bam-lengthfilter:nextflow" : null)
    tag "$meta.RG"
    label 'local'

    input:
    tuple meta, "${meta.RG}.bam" from filterlength_in

    output:
    tuple meta, "${meta.RG}.output.bam" into filterlength_out

    script:
    """
    bam-lengthfilter -c $params.bamfilter_length_cutoff -l $params.compression_level -o \"${meta.RG}.output.bam\" \"${meta.RG}.bam\"
    """
}

process filterLengthCount {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    tag "$meta.RG"
    label 'local'

    input:
    tuple meta, "${meta.RG}.bam" from filterlength_out

    output:
    tuple meta, "${meta.RG}.bam", stdout into filterlengthcount_out
    
    when:
    params.rerun == false   

    script:
    """
    samtools view -c \"${meta.RG}.bam\"
    """
}

//add filterlength count to meta
filterlengthcount_out.map{[it[0]+['ReadsLengthfiltered':it[2].trim()], it[1]]}.into{ gathertaxon_in; tofasta_in }

process toFasta {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    tag "$meta.RG"
    label "local"

    input:
    set meta, "${meta.RG}.bam" from tofasta_in

    output:
    set meta, "${meta.RG}.fa" into tofasta_out

    when:
    meta.ReadsLengthfiltered.toInteger() > 0 && params.rerun == false   

    script:
    """
    samtools fasta \"${meta.RG}.bam\" > \"${meta.RG}.fa\"
    """
}

if (params.testrun){
    // If using test-data, extract the test-database
    process extractTestDatabase {
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"
        label 'process_medium'
        label 'local'
        tag "DB: TestDB"

        input:
        path 'database' from database

        output:
        path 'TestDB' into database_extracted

        when:
        params.rerun == false   

        script:
        """
        tar -xvzf database/database.tar.gz
        """
    }
}

database_out = params.testrun ? database_extracted : database 
tofasta_out.combine(database_out).set{runkraken_in}

process runKrakenUniq {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakenuniq:0.7.3--pl5321h19e8d03_0' :
        'quay.io/biocontainers/krakenuniq:0.7.3--pl5321h19e8d03_0' }"

    publishDir 'stats', mode: 'copy', pattern:"*translate", saveAs: {"${meta.RG}.kraken.translate"}
    publishDir 'stats', mode: 'copy', pattern:"*report", saveAs: {"${meta.RG}.kraken.report"}
    label 'process_high'
    label 'local'
    tag "$meta.RG"

    input:
    tuple meta, "${meta.RG}.fa", "database" from runkraken_in

    output:
    tuple meta, "krakenUniq.translate", "krakenUniq.report" into runkraken_out   

    when:
    params.rerun == false   

    script:
    """
    krakenuniq --threads ${task.cpus} --db database --fasta-input \"${meta.RG}.fa\" --report-file krakenUniq.report > output.krakenUniq
    krakenuniq-translate --db database --mpa-format output.krakenUniq > krakenUniq.translate
    echo "\$(grep -E '^(#|%|[0-9]).*' krakenUniq.report)" > krakenUniq.report
    """
}

process findBestNode{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"
    label 'process_low'
    label 'local'
    tag "$meta.RG"

    input:
    tuple meta, "krakenUniq.translate", "krakenUniq.report" from runkraken_out

    output:
    tuple meta, "parsed_record.tsv", "krakenUniq.translate" into bestspecies_out

    when:
    params.rerun == false   

    script:
    """
    parse_report.py krakenUniq.report ${params.krakenuniq_min_kmers} ${params.krakenuniq_min_reads}
    """
}

taxid.splitCsv(sep:'\t')
    .map{[it[0], it[2]]} //[tax_id, species_file]}
    .unique()
    .groupTuple()
    .set{taxid}
    
bestspecies_out.map{meta, record, translate -> [meta, translate, record.splitCsv(sep:'\t', header:true)]}
    .transpose()
    .map{meta, translate, record -> [record.BestTaxID, meta+record, translate]}
    .combine(taxid,by:0)
    .map{taxid,meta,translate,species -> [meta,translate,species]}
    .set{bestspecies_out}

// This block replaces the default mappings assigned by kraken by those
// specified in the fixed file

def famList = []
specmap = new File("${params.fixed}")
specs = specmap.exists() ? Channel.fromPath("${params.fixed}") : Channel.empty()

if(specmap.exists()){
    specmap.eachLine{famList << it.split("\t").flatten()[0]}
    specs.splitCsv(sep:'\t', header:['fam','sp_tag','path'], skip:1)
         .map{row -> [row.fam, row.sp_tag, file(row.path)]}
         .into{specs; rerun}
} else {
    rerun = Channel.empty()
}    


// if rerun --> get extracted bam files from each family + the report file to replace meta variable
// this channel replaces the mapbwa_in channel downstream, so we need:
// [meta, extracted, ref]

rerun.combine(outdir, by:0)
  .map{fam,sp,ref,bam,rg -> [rg,fam,sp,ref,bam]}
  .combine(report, by:[0,1])
  .map{rg,fam,sp,ref,bam,meta -> [meta+['Reference':'fixed','Species':sp,'ref_path':ref, 'Taxon':meta.ExtractLVL=='f' ? meta.Family : meta.Order], bam, ref]} 
  .unique{it[0].RG+it[0].Family}
  .set{rerun}

bestspecies_out
   .map{meta,translate,species -> [meta.Family,meta,translate,species]}
   .branch {
       replace: it[0] in famList
       keep: true
    }
   .set{bestspecies_out}

bestspecies_out.replace
    .combine(specs, by:0)
    .map{fam,meta,translate,species,sp_tag,path -> [fam, meta+['Reference':'fixed', 'ref_path':path], translate, sp_tag]}
    .set{bestspecies_replaced}

bestspecies_out.keep
    .mix(bestspecies_replaced)
    .map{fam,meta,translate,specs -> [meta.RG, meta, translate, specs]}
    .set{bestspecies_out} 

gathertaxon_in.map{meta, bam -> [meta.RG, meta, bam]}
    .cross(bestspecies_out)
    .map{a,b -> [b[1], a[2], b[2], b[3]]} // [meta, bam, translate, [species]]
    .map{meta,bam,translate,species -> [meta+['Taxon': meta.ExtractLVL == 'f' ? meta.Family : meta.Order ],bam,translate,species]}
    .set{gathertaxon_in}

process gatherByTaxon {
    label 'process_low'
    label 'local'
    tag "$meta.RG:$meta.Taxon"

    input:
    tuple meta, "${meta.RG}.bam", "${meta.RG}.translate", references from gathertaxon_in

    output:
    tuple meta, "${meta.RG}.bam", 'ids.txt', references, stdout into gathertaxon_out

    script:
    """
    grep "${meta.ExtractLVL}__${meta.Taxon}" \"${meta.RG}.translate\" | cut -f1 | tee ids.txt | wc -l
    """
}

gathertaxon_out
    .map{meta,bam,ids,ref,count -> [meta+['ReadsExtracted':count.trim() as int], bam,ids,ref]}
    .set{extractbam_in}

process extractBam {
    container (workflow.containerEngine ? "merszym/bamfilter:nextflow" : null)
    tag "$meta.RG:$meta.Taxon"
    label "process_low"
    label 'local'

    input:
    set meta, "${meta.RG}.bam", 'ids.txt', references from extractbam_in

    output:
    set meta, "${meta.RG}.output.bam", references into extractbam_out

    when:
    meta.ReadsExtracted >= params.krakenuniq_min_reads 

    script:
    """
    bamfilter -i ids.txt -l $params.compression_level -o \"${meta.RG}.output.bam\" \"${meta.RG}.bam\"
    """
}

extractbam_out
    .transpose()
    .map{meta,bam,ref -> [meta+["Species":ref],bam]}
    .set{sortbam_in}

process sortBam{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    publishDir 'out', mode: 'copy', saveAs: {out_bam}
    tag "${meta.RG}:${meta.Taxon}"
    label "process_low"
    label "local"

    input:
    tuple meta, "${meta.Taxon}.extracted.bam" from sortbam_in

    output:
    tuple meta, "${meta.Taxon}.sorted.bam" into sortbam_out

    script:
    out_bam = "${meta.Taxon}/1-extracted/${meta.RG}_extractedReads-${meta.Taxon}.bam"
    """
    samtools sort -n -l $params.compression_level -o \"${meta.Taxon}.sorted.bam\"  \"${meta.Taxon}.extracted.bam\" 
    """
}

sortbam_out
    .combine(genomesdir)
    .map{ meta,bam,genomes -> [meta, bam, meta.ref_path ?: genomes] }
    .set{mapbwa_default}

mapbwa_in = params.rerun ? rerun : mapbwa_default

process mapBwa {
    container (workflow.containerEngine ? "merszym/network-aware-bwa:v0.5.10" : null)
    publishDir 'out', mode: 'copy', saveAs: {out_bam}, pattern: '*.bam'
    tag "${meta.RG}:${meta.Taxon}:${meta.Species}"
    label "process_low"
    label 'local'

    input:
    tuple meta, "${meta.Taxon}.sorted.bam", "reference" from mapbwa_in

    output:
    tuple meta, "${meta.Taxon}.mapped.bam" into mapbwa_out

    script:
    index = meta.Reference=='fixed' ? "bwa index reference" : ""
    genome = meta.Reference=='fixed' ? "reference" : "reference/${meta.Family}/${meta.Species}.fasta"
    out_bam = "${meta.Taxon}/${meta.Reference}/2-aligned/${meta.RG}.${meta.Family}.${meta.Species}.bam"

    """
    $index
    bwa bam2bam -g \"${genome}\" -n 0.01 -o 2 -l 16500 --only-aligned \"${meta.Taxon}.sorted.bam\" > \"${meta.Taxon}.mapped.bam\"
    """
}

process filterMappedBam{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    tag "${meta.RG}:${meta.Taxon}:${meta.Species}"
    label "process_low"
    label "local"

    input:
    tuple meta, "${meta.Taxon}.mapped.bam" from mapbwa_out

    output:
    tuple meta, "${meta.Taxon}.mapped_filtered.bam", stdout into filtermappedbam_out

    script:
    """
    samtools view -b -u -q $params.bamfilter_quality_cutoff \"${meta.Taxon}.mapped.bam\" \
    | samtools sort -l $params.compression_level -o \"${meta.Taxon}.mapped_filtered.bam\"
    samtools view -c \"${meta.Taxon}.mapped_filtered.bam\"
    """
}

filtermappedbam_out
    .map{meta,bam,count -> [meta+['ReadsMapped':count.trim() as int], bam]}
    .set{dedupbam_in}

process dedupBam {
    container (workflow.containerEngine ? "merszym/biohazard_bamrmdup:v0.2.2" : null)
    publishDir 'out', mode: 'copy', pattern: "*.bam", saveAs: {out_bam}
    tag "${meta.RG}:${meta.Family}:${meta.Species}"
    label "process_low"
    label 'local'

    input:
    tuple meta, "${meta.Species}.bam" from dedupbam_in

    output:
    tuple meta, "${meta.Species}.deduped.bam" into dedupedbam_out

    script:
    out_bam = "${meta.Taxon}/${meta.Reference}/3-deduped/${meta.RG}.${meta.Family}.${meta.Species}_deduped.bam"
    """
    bam-rmdup -r -o \"${meta.Species}.deduped.bam\" \"${meta.Species}.bam\" > rmdup.txt
    """
}

process getDedupStats{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    tag "${meta.RG}:${meta.Family}:${meta.Species}"
    label "process_low"
    label "local"

    input:
    tuple meta, "${meta.Species}.deduped.bam" from dedupedbam_out

    output:
    tuple meta, "${meta.Species}.deduped.bam", stdout, "count.txt" into dedupedstats_out

    script:
    """
    samtools coverage -H \"${meta.Species}.deduped.bam\" | cut -f 5
    samtools view -c \"${meta.Species}.deduped.bam\" > count.txt
    """
}

//sort by readgroup, family and if same (?:) by count. GroupTuple by readgroup and family
//continue with the best species per family
dedupedstats_out
    .map{meta,bam,coverage,count -> [meta+['CoveredBP':coverage.trim() as int, 'ReadsDeduped':count.text.strip() as int], bam]}
    .map{meta,bam -> [meta.RG, meta.Family, meta.CoveredBP, meta, bam]}
    .toSortedList({ a,b -> a[0]+a[1] <=> b[0]+b[1] ?: a[2] <=> b[2]})
    .flatMap{n -> n[0..-1]}
    .groupTuple(by:[0,1])   //[[rg, fam, [covered_bp < .. < covered_bp][meta,meta,meta],[bam,bam,bam]]
    .map{n -> [n[3][-1], n[4][-1]]} // [meta, bam]
    .branch{
        no_bed: it[0].Family in famList
        bed: true
    }
    .set{dedupedstats_out}

dedupedstats_out.bed
    .combine(bedfiledir)
    .set{runbed_in} 

process runIntersectBed{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3' :
        'quay.io/biocontainers/bedtools:2.30.0--h468198e_3' }"
    label "process_low"
    label 'local'
    tag "${meta.RG}:${meta.Family}:${meta.Species}"
    publishDir 'out', mode: 'copy', saveAs: {out_bam}

    input:
    tuple meta, "${meta.Species}.bam", "masked" from runbed_in

    output:
    tuple meta, "${meta.Species}.masked.bam" into runbed_out
    
    script:
    out_bam = "${meta.Taxon}/${meta.Reference}/4-bedfiltered/${meta.RG}.${meta.Family}.${meta.Species}_deduped_bedfiltered.bam"
    """
    bedtools intersect -a \"${meta.Species}.bam\" -b masked/\"${meta.Species}.masked.bed\" -v > \"${meta.Species}.masked.bam\"
    """
}

process getBedfilteredCounts{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    tag "${meta.RG}:${meta.Family}:${meta.Species}"
    label "process_low"
    label "local"

    input:
    tuple meta, "${meta.Species}.masked.bam" from runbed_out

    output:
    tuple meta, "${meta.Species}.masked.bam", stdout, 'coverage.txt' into bedfiltercounts_out

    """
    samtools view -c \"${meta.Species}.masked.bam\"
    samtools coverage -H \"${meta.Species}.masked.bam\" | cut -f 5 > 'coverage.txt'
    """
}

bedfiltercounts_out
    .map{meta,bam,count,cov -> [meta+['ReadsBedfiltered':count.trim() as int, 'PostBedCoveredBP':cov.text.trim() as int], bam]}
    .set{bedfiltercounts_out}

//calculate the percentage of a family here
bedfiltercounts_out.map{meta, bam -> [meta, bam, meta.ReadsBedfiltered]}.set{bedfiltercounts_out}

dedupedstats_out.no_bed
    .map{meta,bam -> [meta, bam, meta.ReadsDeduped]}
    .mix(bedfiltercounts_out)
    .into{bedfiltercounts_out;total_rg}

total_rg.map{meta,bam,count -> [meta.RG, count]}
    .groupTuple()    // [RG, [count, count, count, ...]]
    .map{rg,count -> [rg, count.sum()]}
    .set{total_rg}

bedfiltercounts_out.map{meta,bam,count -> [meta.RG, meta, bam, count]}
    .combine(total_rg, by:0)
    .map{rg,meta,bam,count,total_rg -> [meta+['FamPercentage':total_rg==0 ? 0: (count*100/total_rg).trunc(2) ], bam]}
    .map{meta,bam -> [meta+['FamPercentage': params.rerun ? '-' : meta.FamPercentage], bam]} // the calculation doesnt work in rerun-context
    .set{bedfiltercounts_out}

bedfiltercounts_out.branch{
    extract: it[0].Family in famList
    only_stats: true
}.set{damageanalysis_in}

process getDeaminationStats{
    container (workflow.containerEngine ? "merszym/bam_deam:nextflow" : null)
    tag "${meta.RG}:${meta.Taxon}:${meta.Species}"
    label 'process_medium'
    label 'local'
    
    input:
    tuple meta, "${meta.Species}.bam" from damageanalysis_in.only_stats
    
    output:
    tuple meta, 'ancient_stats.tsv' into damagestats_out
    
    script:
    """
    bam_deam_stats.py \"${meta.Species}.bam\" only_stats > ancient_stats.tsv
    """
}

process extractDeaminatedReads{
    container (workflow.containerEngine ? "merszym/bam_deam:nextflow" : null)
    publishDir 'out', mode: 'copy', saveAs: { out_bam1 }, pattern:"*1.bam"
    publishDir 'out', mode: 'copy', saveAs: { out_bam2 }, pattern:"*3.bam"
    tag "${meta.RG}:${meta.Taxon}:${meta.Species}"
    label 'process_medium'
    label 'local'
    
    input:
    tuple meta, "${meta.Species}.bam" from damageanalysis_in.extract
    
    output:
    tuple meta, 'output.deaminated1.bam', 'output.deaminated3.bam', "${meta.Species}.bam" into masking_in
    tuple meta, 'ancient_stats.tsv' into damageanalysis_out
    
    script:
    out_bam1 = "${meta.Taxon}/${meta.Reference}/5-deaminated/${meta.RG}.${meta.Family}.${meta.Species}_deduped_deaminated_1term.bam"
    out_bam2 = "${meta.Taxon}/${meta.Reference}/5-deaminated/${meta.RG}.${meta.Family}.${meta.Species}_deduped_deaminated_3term.bam"
    """
    bam_deam_stats.py \"${meta.Species}.bam\" > ancient_stats.tsv
    """
}

process maskDeamination{
    container (workflow.containerEngine ? "merszym/bam_deam:nextflow" : null)
    tag "${meta.RG}:${meta.Taxon}:${meta.Species}"
    label 'process_medium'
    label 'local'
    
    input:
    tuple meta, 'deaminated1.bam', 'deaminated3.bam', "all_reads.bam" from masking_in
    
    output:
    tuple meta, 'deaminated1.masked.bam', 'deaminated3.masked.bam', "all_reads.bam" into mpileup_in
    
    script:
    """
    mask_qual_scores.py deaminated1.bam
    mask_qual_scores.py deaminated3.bam
    """
}


process createMpileups{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    publishDir {out}, mode: 'copy', pattern: '*.tsv'
    tag "${meta.RG}:${meta.Taxon}:${meta.Species}"
    label 'process_medium'
    label 'local'
    
    input:
    tuple meta, 'in.deaminated1.bam', 'in.deaminated3.bam', "all_reads.bam" from mpileup_in
    
    output:
    tuple "${meta.RG}.${meta.Family}.${meta.Species}_all_mpiled.tsv",
          "${meta.RG}.${meta.Family}.${meta.Species}_term1_mpiled.tsv", 
          "${meta.RG}.${meta.Family}.${meta.Species}_term3_mpiled.tsv" 
   
    script:
    out = "out/${meta.Taxon}/${meta.Reference}/6-mpileups/"
    args = "--output-BP-5 --no-output-ends --no-output-ins --no-output-del"
    """
    samtools mpileup all_reads.bam $args  > \"${meta.RG}.${meta.Family}.${meta.Species}_all_mpiled.tsv\"
    samtools mpileup in.deaminated1.bam $args  > \"${meta.RG}.${meta.Family}.${meta.Species}_term1_mpiled.tsv\"
    samtools mpileup in.deaminated3.bam $args > \"${meta.RG}.${meta.Family}.${meta.Species}_term3_mpiled.tsv\" 
    """
}

damageanalysis_out.mix(damagestats_out)
    .map{meta,damage -> [meta,damage.splitCsv(sep:'\t', header:true)]}
    .transpose()
    .map{meta,dmg -> meta+dmg}
    .set{combine_meta}

final_report = params.rerun ? final_report : Channel.empty()
final_report.map{it[2]}.set{final_report}

//
//
// Write Reports
//
//

combine_meta
  .map{ meta -> meta+[
          'ProportionMapped': ((meta.ReadsExtracted as int)==0 || meta.ReadsMapped == 0) ? 0 : (meta.ReadsMapped/(meta.ReadsExtracted as int)).trunc(4),
          'DuplicationRate': (meta.ReadsDeduped==0 || meta.ReadsMapped==0) ? 0 : (meta.ReadsMapped/meta.ReadsDeduped).trunc(4)
        ]
      }
  .mix(final_report)
  .unique{meta -> meta.RG+meta.Species+meta.Reference}
  .into{
    summary_file ;
    damageanalysis_file;
    bedfiltercounts_file;
    dedupedstats_file;
    filtermappedbam_file;
    gathertaxon_file;
    splitcount_file
  }

header_map = [
 'tax'    : 'Order\tFamily\tSpecies\tReference',
 'split'  : 'ReadsRaw\tReadsFiltered\tReadsLengthfiltered',
 'kraken' : 'FamilyKmers\tKmerCoverage\tKmerDupRate',
 'deam'   : 'Ancientness\tReadsDeam(1term)\tReadsDeam(3term)\tDeam5(95ci)\tDeam3(95ci)\tDeam5Cond(95ci)\tDeam3Cond(95ci)',
 'extract': 'ExtractLVL\tReadsExtracted',
 'map'    : 'ReadsMapped\tProportionMapped',
 'dedup'  : 'ReadsDeduped\tDuplicationRate\tCoveredBP', 
 'bed'    : 'ReadsBedfiltered\tPostBedCoveredBP',  
]

def getVals = {String header, meta, res=[] ->
    header.split('\t').each{res << meta[it]}
    res.join('\t')
}

summary_file
  .collectFile( name:'final_report.tsv', 
    seed:[
      'RG',
      header_map['split'],
      header_map['kraken'], 
      header_map['extract'],
      header_map['tax'], 
      header_map['map'],
      header_map['dedup'], 
      header_map['bed'],
      'FamPercentage', 
      header_map['deam']
    ].join('\t'), storeDir:'.', newLine:true, sort:true
  ){[
      it.RG,
      getVals(header_map['split'],   it),
      getVals(header_map['kraken'],  it),
      getVals(header_map['extract'], it),
      getVals(header_map['tax'],     it),
      getVals(header_map['map'],     it),
      getVals(header_map['dedup'],   it),
      getVals(header_map['bed'],     it),
      it.FamPercentage,
      getVals(header_map['deam'],    it),
    ].join('\t') 
  }
  .subscribe {
    println get_info_msg("Summary reports saved")
  }

damageanalysis_file
  .collectFile( 
     storeDir: 'stats', newLine:true, 
     seed:[
       header_map['tax'], 
       header_map['deam']
     ].join('\t')
  ){[
      "${it.RG}_04_deamination.tsv", [
         getVals(header_map['tax'],  it),
         getVals(header_map['deam'], it)
      ].join('\t')
    ]
  }

bedfiltercounts_file
  .collectFile(
    storeDir: 'stats', newLine:true, 
    seed: [
      header_map['tax'],
      header_map['bed']
    ].join('\t')
  ){[
      "${it.RG}_03_bedfiltered.tsv", 
         [ 
           getVals(header_map['tax'], it),
           getVals(header_map['bed'], it)
         ].join('\t')
    ]
   }

dedupedstats_file
  .collectFile(
    storeDir: 'stats', newLine:true,
    seed: [
      header_map['tax'],
      header_map['dedup']
    ].join('\t')        
  ){[ 
      "${it.RG}_02_deduped.tsv", 
        [
          getVals(header_map['tax'],   it),
          getVals(header_map['dedup'], it)
        ].join('\t')
    ]
   }

filtermappedbam_file
  .collectFile(
    storeDir: 'stats', newLine:true,
    seed: [
      header_map['tax'],
      header_map['map']
    ].join('\t')        
  ){[ 
      "${it.RG}_01_mapped.tsv", 
        [
          getVals(header_map['tax'], it),
          getVals(header_map['map'], it)
        ].join('\t')
    ]
   }

gathertaxon_file
  .unique{meta -> meta.RG+meta.Taxon}
  .collectFile(storeDir: 'stats', seed:'Taxon\tReadsExtracted', newLine:true) {meta ->
    [ "${meta.RG}_00_extracted.tsv", "${meta.Taxon}\t${meta.ReadsExtracted}"]
  }

splitcount_file
  .unique{it.RG}
  .collectFile(
    storeDir: 'stats', newLine:true,
    seed: [
      'RG',
      header_map['split'],
    ].join('\t')        
  ){[ 
      "splitcounts.tsv", 
        [
          it.RG,
          getVals(header_map['split'], it),
        ].join('\t')
    ]
   }

