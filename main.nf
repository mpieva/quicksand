#!/usr/bin/env nextflow

nextflow.enable.dsl = 1

red = "\033[0;31m"
white = "\033[0m"
cyan = "\033[0;36m"
yellow = "\033[0;33m"

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
def add_to_dict(meta, key, val){
    meta[key] = val
    return meta
}
def get_taxon(meta){
    taxon = params.taxlvl == 'f' ? meta.Family : meta.Order
    return taxon
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
inbam      = params.bam      ? Channel.fromPath("${params.bam}",                 checkIfExists:true)  : Channel.empty()
indexfile  = params.rg       ? Channel.fromPath("${params.rg}",                  checkIfExists:true)  : Channel.empty()
splitdir   = params.split    ? Channel.fromPath("${params.split}/*",             checkIfExists:true)  : Channel.empty()
genomesdir = params.genomes  ? Channel.fromPath("${params.genomes}", type:'dir', checkIfExists:true)  : exit_missing_required("--genomes")
database   = params.db       ? Channel.fromPath("${params.db}", type:'dir',      checkIfExists:true)  : exit_missing_required("--db")
bedfiledir = params.bedfiles ? Channel.fromPath("${params.bedfiles}",type:'dir', checkIfExists:true)  : exit_missing_required("--bedfiles") 
specmap    = params.specmap  ? Channel.fromPath("${params.specmap}",             checkIfExists:true)  : Channel.empty()


// Additional File Channels
nobed_families = params.skip_bed.split(',')
taxid = new File("${params.genomes}/taxid_map.tsv").exists() ? Channel.fromPath("${params.genomes}/taxid_map.tsv", type:'file') : Channel.fromPath("${baseDir}/assets/taxid_map_example.tsv", type:'file') 

if (! new File("${params.genomes}/taxid_map.tsv").exists()) {
    log.info get_warn_msg("The file 'taxid_map.tsv' is missing in your genomes dir! Using fallback ${baseDir}/assets/taxid_map_example.tsv")
}    

splitfile = new File("${params.splitfile}").exists() ? Channel.fromPath("${params.splitfile}", type: "file").splitCsv(skip:1) : Channel.empty()


//
//
// Validation
//
//

// Validate user input
if(params.taxlvl !in ['f','o']){
    exit_with_error_msg("ArgumentError","taxlvl must be one of [o, f] not ${params.taxlvl}")
}

if(params.split && (params.bam || params.rg)){
    log.info get_info_msg("Use: nextflow run mpieva/quicksand {--rg FILE --bam FILE | --split DIR}")
    exit_with_error_msg("ArgumentError", "Too many arguments")
} 
if(!params.split && !(params.bam && params.rg)){
    log.info get_info_msg("Use: nextflow run mpieva/quicksand {--rg FILE --bam FILE | --split DIR}")
    exit_with_error_msg("ArgumentError", "Too few arguments")
}

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
    params.bam && params.rg   

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
   .map{[['id': it.baseName], it]}      
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
    tuple meta, "${meta.id}.CC.txt" into crosscont_out 

    script:
    """
    cross_cont.py splittingstats.txt > \"${meta.id}.CC.txt\"
    """
}


crosscont_out
    .filter{it[1].text != ''}
    .collectFile(storeDir: '.') {meta, file ->
        [ "${meta.id}_cc.txt", file.text]
    }


process fastq2Bam{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    tag "$meta.id"
    label 'local'
    
    input:
    tuple meta, "${meta.id}.fq" from splitbam_out.fastq
    
    output:
    tuple meta, "${meta.id}.bam" into fastq2bam_out
    
    script:
    """
    samtools import -0 \"${meta.id}.fq\" -o \"${meta.id}.bam\"
    """
}

splitbam_out.bam.mix(fastq2bam_out).set{filterbam_in}

process filterBam {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    label 'local'
    tag "$meta.id:Flag:${params.bamfilterflag}"

    input:
    tuple meta, "${meta.id}.bam" from filterbam_in

    output:
    tuple meta, "${meta.id}.filtered.bam", stdout, 'filtercount.txt' into filterbam_out

    script:
    """
    samtools view -c -F 128 \"${meta.id}.bam\"
    samtools view -b -u -F ${params.bamfilterflag} -o \"${meta.id}.filtered.bam\" \"${meta.id}.bam\"
    samtools view -c \"${meta.id}.filtered.bam\" > 'filtercount.txt'
    """
}

//add splitcounts to meta
filterbam_out
    .map{[add_to_dict(it[0],'splitcount',it[2].trim()), it[1], it[3].text.trim()]}
    .map{[add_to_dict(it[0],'filtercount',it[2]), it[1]]}
    .set{filterlength_in}

process filterLength {
    container (workflow.containerEngine ? "merszym/bam-lengthfilter:nextflow" : null)
    tag "$meta.id"
    label 'local'

    input:
    tuple meta, "${meta.id}.bam" from filterlength_in

    output:
    tuple meta, "${meta.id}.output.bam" into filterlength_out

    script:
    """
    bam-lengthfilter -c $params.bamfilter_length_cutoff -l $params.compression_level -o \"${meta.id}.output.bam\" \"${meta.id}.bam\"
    """
}

process filterLengthCount {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    tag "$meta.id"
    label 'local'

    input:
    tuple meta, "${meta.id}.bam" from filterlength_out

    output:
    tuple meta, "${meta.id}.bam", stdout into filterlengthcount_out
    
    script:
    """
    samtools view -c \"${meta.id}.bam\"
    """
}

//add filterlength count to meta
filterlengthcount_out.map{[add_to_dict(it[0],'lengthfiltercount',it[2].trim()), it[1]]}.into{ gathertaxon_in; splitcount_file; tofasta_in }

// Include the old splitcounts file in the summary
splitcount_file.map{meta,bam -> ["$meta.id\t$meta.splitcount\t$meta.filtercount\t$meta.lengthfiltercount"]}
    .concat(splitfile)
    .unique{it[0]}
    .map{it[0]}
    .collectFile(storeDir: 'stats', name: "splitcounts.tsv", newLine: true,
                 seed: "RG\tReadsRaw\tReadsFiltered\tReadsLengthfiltered")

process toFasta {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    tag "$meta.id"
    label "local"

    input:
    set meta, "${meta.id}.bam" from tofasta_in

    output:
    set meta, "${meta.id}.fa" into tofasta_out

    when:
    meta.lengthfiltercount.toInteger() > 0

    script:
    """
    samtools fasta \"${meta.id}.bam\" > \"${meta.id}.fa\"
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

    publishDir 'stats', mode: 'copy', pattern:"*translate", saveAs: {"${meta.id}.kraken.translate"}
    publishDir 'stats', mode: 'copy', pattern:"*report", saveAs: {"${meta.id}.kraken.report"}
    label 'process_high'
    label 'local'
    tag "$meta.id"

    input:
    tuple meta, "${meta.id}.fa", "database" from runkraken_in

    output:
    tuple meta, "krakenUniq.translate", "krakenUniq.report" into runkraken_out   

    script:
    """
    krakenuniq --threads ${task.cpus} --db database --fasta-input \"${meta.id}.fa\" --report-file krakenUniq.report > output.krakenUniq
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
    tag "$meta.id"

    input:
    tuple meta, "krakenUniq.translate", "krakenUniq.report" from runkraken_out

    output:
    tuple meta, "parsed_record.tsv", "krakenUniq.translate" into bestspecies_out

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
// specified in the specmap file

def famList = []
specmap = new File("${params.specmap}")
specs = specmap.exists() ? Channel.fromPath("${params.specmap}") : Channel.empty()

if(specmap.exists()){
    specmap.eachLine{famList << it.split("\t").flatten()[0]}
    specs.splitCsv(sep:'\t')
         .map{fam, sp -> [fam, sp.split(',').flatten()]}
         .set{specs}
}     

bestspecies_out
   .map{meta,translate,species -> [meta.Family,meta,translate,species]}
   .branch {
       replace: it[0] in famList
       keep: true
    }
   .set{bestspecies_out}

bestspecies_out.replace
    .combine(specs, by:0) //[family, meta, translate, species, [new_species,...])]
    .map{fam,meta,translate,species,new_species -> [fam, meta, translate, new_species]}
    .set{bestspecies_replaced}

bestspecies_out.keep
    .mix(bestspecies_replaced)
    .map{fam,meta,translate,specs -> [meta.id, meta, translate, specs]}
    .set{bestspecies_out} 

gathertaxon_in.map{meta, bam -> [meta.id, meta, bam]}
    .cross(bestspecies_out)
    .map{a,b -> [b[1], a[2], b[2], b[3]]} // [meta, bam, translate, [species]]
    .map{meta,bam,translate,species -> [add_to_dict(meta,'Taxon',get_taxon(meta)),bam,translate,species]}
    .set{gathertaxon_in}

process gatherByTaxon {
    label 'process_low'
    label 'local'
    tag "$meta.id:$meta.Taxon"

    input:
    tuple meta, "${meta.id}.bam", "${meta.id}.translate", references from gathertaxon_in

    output:
    tuple meta, "${meta.id}.bam", 'ids.txt', references, stdout into (gathertaxon_out, gathertaxon_file)

    script:
    """
    grep "${params.taxlvl}__${meta.Taxon}" \"${meta.id}.translate\" | cut -f1 | tee ids.txt | wc -l
    """
}

gathertaxon_out
    .map{meta,bam,ids,ref,count -> [add_to_dict(meta,'Extracted',count.trim() as int), bam,ids,ref]}
    .set{extractbam_in}

gathertaxon_file
        .unique{it[0].id+it[0].Taxon}
        .collectFile(storeDir: 'stats', seed:'Taxon\tReadsExtracted', newLine:true) {meta, bam, ids, ref, count ->
            [ "${meta.id}_extracted.tsv", "${meta.Taxon}\t${count.trim()}"]
        }

process extractBam {
    container (workflow.containerEngine ? "merszym/bamfilter:nextflow" : null)
    publishDir 'out', mode: 'copy', saveAs: {out_bam}
    tag "$meta.id:$meta.Taxon"
    label "process_low"
    label 'local'

    input:
    set meta, "${meta.id}.bam", 'ids.txt', references from extractbam_in

    output:
    set meta, "${meta.id}.output.bam", references into extractbam_out

    when:
    meta.Extracted >= params.krakenuniq_min_reads 

    script:
    out_bam = params.byrg ? "${meta.id}/${meta.Taxon}/unmapped/${meta.id}_extractedReads-${meta.Taxon}.bam" : "${meta.Taxon}/unmapped/${meta.id}_extractedReads-${meta.Taxon}.bam"
    
    """
    bamfilter -i ids.txt -l $params.compression_level -o \"${meta.id}.output.bam\" \"${meta.id}.bam\"
    """
}

extractbam_out
    .transpose()
    .map{meta,bam,ref -> [add_to_dict(meta,"Species",ref),bam]}
    .set{sortbam_in}

process sortBam{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    tag "${meta.id}:${meta.Taxon}"
    label "process_low"
    label "local"

    input:
    tuple meta, "${meta.Taxon}.extracted.bam" from sortbam_in

    output:
    tuple meta, "${meta.Taxon}.sorted.bam" into sortbam_out

    script:
    """
    samtools sort -n -l $params.compression_level -o \"${meta.Taxon}.sorted.bam\"  \"${meta.Taxon}.extracted.bam\" 
    """
}

sortbam_out
    .combine(genomesdir)
    .set{mapbwa_in}

process mapBwa {
    container (workflow.containerEngine ? "merszym/network-aware-bwa:v0.5.10" : null)
    publishDir 'out', mode: 'copy', saveAs: {out_bam}, pattern: '*.bam'
    tag "${meta.id}:${meta.Taxon}:${meta.Species}"
    label "process_low"
    label 'local'

    input:
    tuple meta, "${meta.Taxon}.sorted.bam", "genomes" from mapbwa_in

    output:
    tuple meta, "${meta.Taxon}.mapped.bam" into mapbwa_out

    script:
    out_bam = params.byrg ? "${meta.id}/${meta.Taxon}/aligned/${meta.Family}.${meta.Species}.bam" : "${meta.Taxon}/aligned/${meta.id}.${meta.Family}.${meta.Species}.bam"

    """
    bwa bam2bam -g genomes/\"${meta.Family}\"/\"${meta.Species}.fasta\" -n 0.01 -o 2 -l 16500 --only-aligned \"${meta.Taxon}.sorted.bam\" > \"${meta.Taxon}.mapped.bam\"
    """
}

process filterMappedBam{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    tag "${meta.id}:${meta.Taxon}:${meta.Species}"
    label "process_low"
    label "local"

    input:
    tuple meta, "${meta.Taxon}.mapped.bam" from mapbwa_out

    output:
    tuple meta, "${meta.Taxon}.mapped_filtered.bam", stdout into (filtermappedbam_out, filtermappedbam_file)

    script:
    """
    samtools view -b -u -q $params.bamfilter_quality_cutoff \"${meta.Taxon}.mapped.bam\" \
    | samtools sort -l $params.compression_level -o \"${meta.Taxon}.mapped_filtered.bam\"
    samtools view -c \"${meta.Taxon}.mapped_filtered.bam\"
    """
}

filtermappedbam_file
    .collectFile(storeDir: 'stats', seed:'Order\tFamily\tSpecies\tReadsMapped', newLine:true) {meta,bam,count ->
        [ "${meta.id}_mapped.tsv", "${meta.Order}\t${meta.Family}\t${meta.Species}\t${count.trim()}"]
    }
filtermappedbam_out
    .map{meta,bam,count -> [add_to_dict(meta,'Mapped',count.trim() as int), bam]}
    .set{dedupbam_in}

process dedupBam {
    container (workflow.containerEngine ? "merszym/biohazard_bamrmdup:v0.2.2" : null)
    publishDir 'out', mode: 'copy', pattern: "*.bam", saveAs: {out_bam}
    tag "${meta.id}:${meta.Family}:${meta.Species}"
    label "process_low"
    label 'local'

    input:
    tuple meta, "${meta.Species}.bam" from dedupbam_in

    output:
    tuple meta, "${meta.Species}.deduped.bam" into dedupedbam_out

    script:
    out_bam = params.byrg ? "${meta.id}/${meta.Taxon}/deduped/${meta.Family}.${meta.Species}_deduped.bam" : "${meta.Taxon}/deduped/${meta.id}.${meta.Family}.${meta.Species}_deduped.bam"
    """
    bam-rmdup -r -o \"${meta.Species}.deduped.bam\" \"${meta.Species}.bam\" > rmdup.txt
    """
}

process getDedupStats{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    tag "${meta.id}:${meta.Family}:${meta.Species}"
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

dedupedstats_out
    .map{meta,bam,coverage,count -> [add_to_dict(meta,'CoveredBP',coverage.trim() as int), bam, count]}
    .map{meta,bam,count -> [add_to_dict(meta,'Deduped', count.text.strip() as int), bam]}
    .into{dedupedstats_out;dedupedstats_file;dedupedstats_file_counts}

//write stats to file
dedupedstats_file
    .collectFile(storeDir: 'stats',newLine:true,seed:'Order\tFamily\tSpecies\tCoveredBP') {meta,bam ->
        [ "${meta.id}_mapped_coverage.tsv", "${meta.Order}\t${meta.Family}\t${meta.Species}\t${meta.CoveredBP}"]
    }

dedupedstats_file_counts
    .collectFile(storeDir: 'stats', newLine:true,seed:'Order\tFamily\tSpecies\tReadsDeduped') {meta,bam ->
        [ "${meta.id}_mapped_deduped.tsv", "${meta.Order}\t${meta.Family}\t${meta.Species}\t${meta.Deduped}"]
    }

//continue with pipeline
dedupedstats_out
    //sort by readgroup, family and if same (?:) by count. GroupTuple by readgroup and family
    //continue with the best species per family
    .map{meta,bam -> [meta.id, meta.Family, meta.CoveredBP, meta, bam]}
    .toSortedList({ a,b -> a[0]+a[1] <=> b[0]+b[1] ?: a[2] <=> b[2]})
    .flatMap{n -> n[0..-1]}
    .groupTuple(by:[0,1])   //[[rg, fam, [covered_bp < .. < covered_bp][meta,meta,meta],[bam,bam,bam]]
    .map{n -> [n[3][-1], n[4][-1]]} // [meta, bam]
    .branch{
        no_bed: it[0].Family in nobed_families
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
    tag "${meta.id}:${meta.Family}:${meta.Species}"
    publishDir 'out', mode: 'copy', saveAs: {out_bam}

    input:
    tuple meta, "${meta.Species}.bam", "masked" from runbed_in

    output:
    tuple meta, "${meta.Species}.masked.bam" into runbed_out
    
    script:
    out_bam = params.byrg ? "${meta.id}/${meta.Taxon}/bedfiltered/${meta.Family}.${meta.Species}_deduped_bedfiltered.bam" : 
                            "${meta.Taxon}/bedfiltered/${meta.id}.${meta.Family}.${meta.Species}_deduped_bedfiltered.bam"
    """
    bedtools intersect -a \"${meta.Species}.bam\" -b masked/\"${meta.Species}.masked.bed\" -v > \"${meta.Species}.masked.bam\"
    """
}

process getBedfilteredCounts{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    tag "${meta.id}:${meta.Family}:${meta.Species}"
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
    .map{meta,bam,count,cov -> [add_to_dict(meta, 'Bedfiltered', count.trim() as int), bam, cov]}
    .map{meta,bam,cov -> [add_to_dict(meta, 'PostBedCoverage', cov.text.trim() as int), bam]}
    .into{bedfiltercounts_file;bedfiltercounts_out}

bedfiltercounts_file
    .collectFile(storeDir: 'stats', newLine:true, seed:'Order\tFamily\tSpecies\tReadsBedfiltered\tPostBedCoveredBP') {meta, bam ->
        [ "${meta.id}_mapped_deduped_bedfiltered.tsv", "${meta.Order}\t${meta.Family}\t${meta.Species}\t${meta.Bedfiltered}\t${meta.PostBedCoverage}"]
    }

//calculate the percentage of a family here
bedfiltercounts_out.map{meta, bam -> [meta, bam, meta.Bedfiltered]}.set{bedfiltercounts_out}

dedupedstats_out.no_bed
    .map{meta,bam -> [meta, bam, meta.Deduped]}
    .mix(bedfiltercounts_out)
    .into{bedfiltercounts_out;total_rg}

total_rg.map{meta,bam,count -> [meta.id, count]}
    .groupTuple()    // [RG, [count, count, count, ...]]
    .map{rg,count -> [rg, count.sum()]}
    .set{total_rg}

bedfiltercounts_out.map{meta,bam,count -> [meta.id, meta, bam, count]}
    .combine(total_rg, by:0)
    .map{rg,meta,bam,count,total_rg -> [add_to_dict(meta,'FamPercentage',total_rg==0?0:(count*100/total_rg).trunc(2)), bam]}
    .into{bedfiltercounts_out; analysis_skip}

damageanalysis_in = params.skip_analyze ? Channel.empty() : bedfiltercounts_out

process analyzeDeamination{
    container (workflow.containerEngine ? "merszym/bam_deam:nextflow" : null)
    publishDir 'out', mode: 'copy', saveAs: { out_bam1 }, pattern:"*1.bam"
    publishDir 'out', mode: 'copy', saveAs: { out_bam2 }, pattern:"*3.bam"
    tag "${meta.id}:${meta.Taxon}:${meta.Species}"
    label 'process_medium'
    label 'local'
    
    input:
    tuple meta, "${meta.Species}.bam" from damageanalysis_in
    
    output:
    tuple 'output.deaminated1.bam', 'output.deaminated3.bam'
    tuple meta, 'ancient_stats.tsv' into damageanalysis_out
    
    script:
    out_bam1 = params.byrg ? "${meta.id}/${meta.Taxon}/deaminated/${meta.Family}.${meta.Species}_deduped_deaminated_1term.bam" : 
                            "${meta.Taxon}/deaminated/${meta.id}.${meta.Family}.${meta.Species}_deduped_deaminated_1term.bam"
    out_bam2 = params.byrg ? "${meta.id}/${meta.Taxon}/deaminated/${meta.Family}.${meta.Species}_deduped_deaminated_3term.bam" : 
                            "${meta.Taxon}/deaminated/${meta.id}.${meta.Family}.${meta.Species}_deduped_deaminated_3term.bam"
    """
    bam_deam_stats.py \"${meta.Species}.bam\" > ancient_stats.tsv
    """
}

damageanalysis_out.map{meta,damage -> [meta,damage.splitCsv(sep:'\t', header:true)]}
    .transpose()
    .map{meta,dmg -> meta+dmg}
    .into{damageanalysis_out;damageanalysis_file}

damageanalysis_file
    .collectFile( storeDir: 'stats', newLine:true, seed:[
       'Order','Family','Species','Ancientness','ReadsDeam(1term)',
       'ReadsDeam(3term)','Deam5(95ci)','Deam3(95ci)','Deam5Cond(95ci)',
       'Deam3Cond(95ci)'].join('\t')){ 
       meta -> [
                 "${meta.id}_mapped_deduped_deamination.tsv",
                 [ meta.Order,
                   meta.Family,
                   meta.Species,
                   meta.Ancientness,
                   meta.ReadsDeam1,
                   meta.ReadsDeam3,
                   meta.Deam5P,
                   meta.Deam3P,
                   meta.Deam5Cond,
                   meta.Deam3Cond
                 ].join('\t')
               ]
   }

//now the final report
if(! params.skip_report){
    analysis_skip.map{meta, bam -> meta}.set{analysis_skip}

    report_in = params.skip_analyze ? analysis_skip : damageanalysis_out
    //add vals that are not in meta
    report_in.map{ meta -> [
    	meta,
    	prop_mapped = (meta.Extracted==0 || meta.Mapped == 0) ? 0 : (meta.Mapped/meta.Extracted).trunc(4),
    	dupl_rate = (meta.Deduped==0 || meta.Mapped==0) ? 0 : (meta.Mapped/meta.Deduped).trunc(4)
        ]
    }
    .collectFile(
        seed:[
                'RG','FamilyKmers','KmerCoverage',
                'KmerDupRate','Order','Family',
                'Species','ReadsRaw','ReadsFiltered','ReadsLengthfiltered','ExtractLVL','ReadsExtracted',
                'ReadsMapped','ProportionMapped','ReadsDeduped',
                'DuplicationRate','CoveredBP','ReadsBedfiltered','PostBedCoveredBP',
                'FamPercentage','Ancientness','ReadsDeam(1term)','ReadsDeam(3term)',
                'Deam5(95ci)','Deam3(95ci)','Deam5Cond(95ci)','Deam3Cond(95ci)'
        ].join('\t'),
        storeDir:'.', newLine:true, sort:true
    ){ meta, pm, dup -> [
            'final_report.tsv', [
                meta.id,
                meta.FamKmers,
                meta.FamKmerCov,
                meta.FamKmerDup,
                meta.Order,
                meta.Family,
                meta.Species,
                meta.splitcount,
                meta.filtercount,
                meta.lengthfiltercount,
                params.taxlvl,
                meta.Extracted,
                meta.Mapped,
                pm,
                meta.Deduped,
                dup,
                meta.CoveredBP,
                meta.Bedfiltered ?: '-',
                meta.PostBedCoverage ?: '-',
                meta.FamPercentage,
                meta.Ancientness ?: '-',
                meta.ReadsDeam1 ?: '-',
                meta.ReadsDeam3 ?: '-',
                meta.Deam5P ?: '-',
                meta.Deam3P ?: '-',
                meta.Deam5Cond ?: '-',
                meta.Deam3Cond ?: '-'
            ].join('\t') 
    ]}
}





