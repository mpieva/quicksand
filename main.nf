#!/usr/bin/env nextflow

red = "\033[0;31m"
white = "\033[0m"
cyan = "\033[0;36m"
yellow = "\033[0;33m"

log.info """
[Sediment_nf]: Execution started: ${workflow.start.format('dd.MM.yyyy HH:mm')} ${cyan}

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
// Help
//
//

params.help = false
if (params.help) {
    print file("$baseDir/assets/help.txt").text
    exit 0
}

//
//
// Variable definition and validation
//
//

// Inhouse programs
params.bwa            = '/home/public/usr/bin/bwa'
params.bamrmdup       = '/home/bioinf/usr/bin/bam-rmdup'
params.kraken_uniq    = '/home/merlin_szymanski/bin/Krakenuniq/'

// User-independent params
params.cutoff         = 35
params.quality        = 25
params.min_kmers      = 129    // min 150bp covered (krakenuniq)
params.min_reads      = 3      // min 3 reads assigned (krakenuniq)
params.krakenthreads  = 4      // number of threads per kraken process
params.level          = 0      // bgzf compression level for intermediate files, 0..9
params.keeppaired     = false  // keep paired reads
params.filterunmapped = false  // filter out unmapped (in case of pre-mapping)
params.analyze        = ""
params.report         = ""
params.splitfile      = "stats/splitcounts.tsv"

def env = System.getenv()

def validate_dir(path, flag){
    File test = new File("$path")
    if (!(test.exists() && test.isDirectory())){
        log.info "[quicksand]: ${red}InputError: The $flag directory doesn't exist ${white}\nPath: $path"
        exit 0;
    }
}

def has_ending(file, extension){
    extension.any{
        file.toString().toLowerCase().endsWith(it)
    }
}


// Validate: BAM, RG and SPLIT
params.bam            = ''
params.rg             = ''
params.split          = ''

if(params.split && (params.bam || params.rg)){
    log.info """
    [quicksand]: ${red}ArgumentError: Too many arguments ${white}
    Use: nextflow run mpieva/quicksand {--rg FILE --bam FILE | --split DIR }
    """.stripIndent()
    exit 0
} 
if(!params.split && !(params.bam && params.rg)){
    log.info """
    [quicksand]: ${red}ArgumentError: Too few arguments ${white}
    See: nextflow run mpieva/quicksand --help
    """.stripIndent()
    exit 0
}

inbam = params.bam? Channel.fromPath("${params.bam}") : Channel.empty()
indexfile = params.rg? Channel.fromPath("${params.rg}") : Channel.empty()


// Validate: GENOME and TAXIDMAP
params.genome = env["QS_GENOME"]
if(!params.genome){
    log.info """[quicksand]: ${red}ArgumentError: Missing --genome flag ${white}"""
    exit 0;
} else { 
    validate_dir("${params.genome}","--genome") 
}

if(new File("${params.genome}/taxid_map.tsv").exists()){
    taxid = Channel.fromPath("${params.genome}/taxid_map.tsv", type:'file')
} else {
    taxid = Channel.fromPath("${baseDir}/assets/taxid_map_example.tsv", type:'file')
    log.info """
        [quicksand]: ${yellow}CRITICAL WARNING: The file 'taxid_map.tsv' is missing in your genomes dir!
        In this run you use the fallback ${baseDir}/assets/taxid_map_example.tsv 
        based on the Refseq release 202. This might lead to WRONG OR MISSING assignments. 
        Please check the Docs, how to set up an own taxid_map.tsv ${white}
        """.stripIndent()
}

// Validate: KRAKEN DB
params.db = env["QS_DB"]
if(!params.db){
    log.info """[quicksand]: ${red}ArgumentError: Missing --db flag ${white}"""
    exit 0;
} else {
    validate_dir("${params.db}", "--db")
}


// Validate: BEDFILES
params.bedfiles = env["QS_BEDFILES"]
if(!params.bedfiles){
    log.info """[quicksand]: ${red}ArgumentError: Missing --bedfiles flag ${white}"""
    exit 0;
} else {
    validate_dir("${params.bedfiles}", "--bedfiles")
}

// Optional
params.specmap = env["QS_SPECMAP"]
params.capture = false
params.byrg    = ''

capture_families = []
if(params.capture){
    capture_families = params.capture.split(',')
}

//
//
// Finally: The Pipeline
//
//

process splitBam {
    maxForks 1
    publishDir 'split', mode: 'copy'
    label 'local'

    input:
    file 'input.bam' from inbam
    file 'indices.tsv' from indexfile

    output:
    file '*.bam' into splitfiles mode flatten
    stdout into splitscriptstats

    when:
    params.bam && params.rg   

    script:
    """
    splitbam -s -c $params.level -f indices.tsv --minscore 10 --maxnumber 0 input.bam
    """
}

splitscriptstats
    .collectFile(storeDir: 'stats', name: "splitstats.tsv", newLine: true)

// If split is defined, start the pipeline here
if(params.split){
    File split_test = new File("${params.split}")
    if (!(split_test.exists() && split_test.isDirectory())){
        log.info "[quicksand]: ${red}InputError: The --split directory doesn't exist ${white}\nPath: ${params.split}"
        exit 0;
    }
    Channel.fromPath("${params.split}/*")
        .ifEmpty{error "[quicksand]: InputError: The Split-Directory is empty"}
        .map{ [it.baseName, it] }
        .set{ splitfiles}
} else {
    splitfiles
        .map{ [it.baseName, it] }      
        .set{ splitfiles}
}

//handle fastq-input
splitfiles
    .branch{
        bam: it[1].getExtension() == "bam"
        fastq: has_ending(it[1], ["fastq"])
        fail: true
    }
    .set {splitfiles}

splitfiles.fail
    .view{"[quicksand]: ${yellow}WARNING: ${it[1]} omitted. File has neither bam nor fastq-ending!${white}"}

process fastq2Bam{
    tag "$rg"
    
    input:
    set rg, file(fastqfile) from splitfiles.fastq
    
    output:
    set rg, "output.bam" into converted_bams
    
    script:
    """
    fastq2bam -1 $fastqfile -o output.bam
    """
}

splitfiles.bam.mix(converted_bams)
    .ifEmpty{error "----\n${white}[quicksand]:${red}WorkflowError: No input-files. Scheck SPLIT-dir or RG-BAM combination. Exit pipeline${white}"}
    .into{splitfiles; splitstats}

process splitStats {
    tag "$rg"

    input:
    set rg,'input.bam' from splitstats

    output:
    set rg, stdout into splitcounts

    script:
    """
    samtools view -c input.bam
    """
}

// if keeppaired==True, use an empty channel, else use splitfiles
filter_paired_in = params.keeppaired ? Channel.empty() : splitfiles

process filterPaired {
    tag "$rg"

    input:
    set rg, 'input.bam' from filter_paired_in

    output:
    set rg, 'output.bam' into filter_paired_out

    script:
    """
    samtools view -b -u -F 1 -o output.bam input.bam
    """
}
//here the paths come together again
post_filter_paired = params.keeppaired ? splitfiles.bam : filter_paired_out

// and do the same with the filter-unmapped step
filter_unmapped_in = params.filterunmapped ? post_filter_paired : Channel.empty()

process filterUnmapped {
    tag "$rg"

    input:
    set rg, 'input.bam' from filter_unmapped_in

    output:
    set rg, 'output.bam' into filter_unmapped_out

    script:
    """
    samtools view -b -u -F 4 -o output.bam input.bam
    """
}

post_filter_unmapped = params.filterunmapped ? filter_unmapped_out : post_filter_paired

process filterLength {
    tag "$rg"

    input:
    set rg, 'input.bam' from post_filter_unmapped

    output:
    set rg, 'output.bam', stdout into tofasta_in
    set rg, 'output.bam' into for_extraction
    set rg, stdout into filtercounts

    script:
    """
    bam-lengthfilter -c $params.cutoff -l $params.level -o output.bam input.bam
    samtools view -c output.bam
    """
}

// In case one needs to run several splitted bam-files in batches, dont overwrite
// the old splitcounts file, but gather everything there
if(new File("$params.splitfile").exists()){
    old_splitcount = Channel.fromPath("${params.splitfile}", type: "file")
        .splitCsv(sep:'\t')
        .filter{it[0] != "readgroup"}
} else {
    old_splitcount = Channel.empty()
}

filtercounts.join(splitcounts)
    .concat(old_splitcount)
    .unique{it[0]}
    .map { rg, fc, sc -> "${rg}\t${sc.trim()}\t${fc.trim()}"}
    .collectFile(storeDir: 'stats', name: "splitcounts.tsv", newLine: true,
                 seed: "readgroup\tsplit count\tfiltered count")

process toFasta {
    tag "$rg"

    input:
    set rg, 'input.bam', filtered_count from tofasta_in

    output:
    set rg, 'output.fa' into tofasta_out

    when:
    filtered_count.toInteger() > 0

    script:
    """
    samtools fasta input.bam > output.fa
    """
}

kraken_db = Channel.fromPath("${params.db}", type:"dir")
tofasta_out.combine(kraken_db).set{pre_kraken}

process runKrakenUniq {
	publishDir 'kraken', mode: 'copy', pattern:"*translate", saveAs: {"${rg}.translate"}
    publishDir 'kraken', mode: 'copy', pattern:"*report", saveAs: {"${rg}.report"}
    cpus "${params.krakenthreads}"
    memory '16GB'
    label 'bigmem'
    label 'local'
    tag "$rg"

    input:
    set rg, "input.fa", db from pre_kraken

    output:
    set rg, "krakenUniq.translate" into kraken_assignments
    set rg, "krakenUniq.report" into find_best

    script:
    """
    ${params.kraken_uniq}/krakenuniq --threads ${task.cpus} --db ${db} --fasta-input input.fa --report-file krakenUniq.report > output.krakenUniq
    ${params.kraken_uniq}/krakenuniq-translate --db ${db} --mpa-format output.krakenUniq > krakenUniq.translate
    """
}

process findBestSpecies{
    tag "$rg"

    input:
    set rg, "krakenUniq.report" from find_best

    output:
    set rg, "parsed_record.tsv" into best_species

    script:
    """
    parse_report.py krakenUniq.report ${params.min_kmers} ${params.min_reads}
    """
}

taxid.splitCsv()
    .map{it[0].split('\t').flatten()}
    .map{[it[0], it[2]]}
    .groupTuple()
    .set{taxid}
    
best_species
    .map{[it[0], it[1].readLines()]}
    .transpose()
    .map{[it[1].split('\t')[1], it[0], it[1].split('\t')[0]]}
    .combine(taxid, by:0)  //[taxid, readgroup, family, [species, species]]
    .set{best_species}

// This block replaces the default mappings assigned by kraken by those
// specified in the specmap file

def famList = []
if (new File("${params.specmap}").exists()){
    new File("${params.specmap}").eachLine{famList << it.split("\t").flatten()[0]}

    specs = Channel.fromPath("${params.specmap}", type: "file")
        .splitCsv(sep:'\t')
        .map{[it[1].split(","), it[0]]} //[[species,species], Family]
   
    best_species
        .map{[it[1]+it[2], it[2], it[3]]} //[newKey, family, [species, species]]
        .branch {
            replace: it[1] in famList
            keep: true
        }
        .set{best_species}

    best_species.replace
        .combine(specs, by:1) //[family, newKey, [original_species, ...],[new_species,...])]
        .map{[it[1], it[0], it[3].flatten()]} //[newKey, family, [new_species,...]]
        .set{replace}

    best_species.keep.mix(replace)
        .map{[it[0], it[2]]} // [newKey, [species, ...]]
        .transpose() // [newKey, species]
        .set{best_species_post}

} else {
    best_species
        .transpose()
        .map{[it[1]+it[2], it[3]]}
        .set{best_species_post}
}

for_extraction
    .cross(kraken_assignments)
    .map { [it[0][0], it[0][1], it[1][1], it[1][1].readLines()] }
    .branch {
        assigned_taxa: it[3].any{ it =~ /c__Mammalia.*f__./  }
        empty: true
    }
    .set{for_extraction}

for_extraction.empty
    .view{"[quicksand]: ${yellow}Info: No Kraken-assignments for Readgroup ${it[0]}${white}"}
    .collectFile(storeDir: 'stats') { it -> [ "${it[0]}_extracted.tsv", "\t"] }

for_extraction.assigned_taxa
    .transpose()
    .filter { it[3] =~ /c__Mammalia.*f__./ }
    .map { rg, bam, kraken, asn -> [rg, bam, kraken, (asn =~ /f__([^|]*)/)[0][1]] }
    .unique()
    .ifEmpty{error "----\n${white}[quicksand]:${red} WorkflowError: No families assigned by Kraken at all. Check Input and Database! Exit pipeline${white}"}
    .set { for_extraction }

process gatherByFamily {
    tag "$rg:$family"

    input:
    set rg, 'input.bam', 'kraken.translate', family from for_extraction

    output:
    set family, rg, 'input.bam', 'ids.txt', stdout into (prepared_for_extraction, count_for_stats, extraction_data)

    script:
    """
    grep "c__Mammalia.*f__$family" kraken.translate | cut -f1 | tee ids.txt | wc -l
    """
}

count_for_stats
        .collectFile(storeDir: 'stats') { family, rg, bamf, idf, count ->
            [ "${rg}_extracted.tsv", "${family}\t${count}"]
        }

process extractBam {
    publishDir 'out', mode: 'copy', saveAs: {out_bam}
    tag "$rg:$family"

    input:
    set family, rg, 'input.bam', 'ids.txt', idcount from prepared_for_extraction

    output:
    set rg, family, 'output.bam' into extracted_reads

    when:
    idcount.toInteger() >= params.min_reads 

    script:
    if(params.byrg){
        out_bam = "${rg}/${rg}_extractedReads-${family}.bam"
    } else {
        out_bam = "${family}/${rg}_extractedReads-${family}.bam"
    }
    """
    bamfilter -i ids.txt -l $params.level -o output.bam input.bam
    """
}

extracted_reads
    .map{[it[0]+it[1], it[0], it[1], it[2]]}
    //extracted reads --> [new_Key, rg, fam, ExtractedReads_Hominidae.bam]
    //best_species --> [New_key, species]
    .cross(best_species_post)
    .map{x,y -> [x[1], x[2], y[1], x[3]] }
    //[Readgroup, Hominidae, Homo_sapiens, ExtractedReads_Hominidae.bam]
    .set{extracted_reads}

reference = Channel.fromPath("${params.genome}", type:"dir")
extracted_reads.combine(reference).set{pre_bwa}

process mapBwa {
    publishDir 'out', mode: 'copy', saveAs: { out_bam }, pattern: '*.bam'
    tag "$rg:$family:$species"

    input:
    set rg, family, species, "input.bam", genomes from pre_bwa

    output:
    set family, rg, species, 'output.bam' into mapped_bam
    set family, rg, species, stdout into (mapped_count, mapping_data)

    script:
    if(params.byrg){
        out_bam = "${rg}/aligned/${family}.${species}.bam"
    } else {
        out_bam = "${family}/aligned/${rg}.${species}.bam"
    }
    """
    samtools sort -n -l0 input.bam \
    | $params.bwa bam2bam -g ${genomes}/$family/\"${species}.fasta\"  -n 0.01 -o 2 -l 16500 --only-aligned - \
    | samtools view -b -u -q $params.quality \
    | samtools sort -l $params.level -o output.bam
    samtools view -c output.bam
    """
}

mapped_count
        .collectFile(storeDir: 'stats') { family, rg, species, count ->
            [ "${rg}_mapped.tsv", "${family}\t${species}\t${count}"]
        }


process dedupBam {
    publishDir 'out', mode: 'copy', pattern: "*.bam", saveAs: {out_bam}
    tag "$rg:$family:$species"

    input:
    set family, rg, species, 'input.bam' from mapped_bam

    output:
    set species, family, rg, "output.bam", stdout, "count.txt" into deduped_bam
    set family, rg, species, stdout into (coverage_count, coverage_data)
    set family, rg, species, "count.txt" into (deduped_count, deduped_data)

    script:
    if(params.byrg){
        out_bam = "${rg}/aligned/${family}.${species}_deduped.bam"
    } else {
        out_bam = "${family}/aligned/${rg}.${species}_deduped.bam"
    }
    """
    $params.bamrmdup -r -o output.bam input.bam > rmdup.txt
    samtools coverage -H output.bam | cut -f 5
    samtools view -c output.bam > count.txt
    """
}

coverage_count
    .collectFile(storeDir: 'stats') { family, rg, species, coverage ->
        [ "${rg}_mapped_coverage.tsv", "${family}\t${species}\t${coverage}"]
    }

deduped_count
    .collectFile(storeDir: 'stats') { family, rg, species, count_file ->
        [ "${rg}_unique_mapped.tsv", "${family}\t${species}\t" + count_file.text]
    }

deduped_bam
    /*sort by readgroup, family and if same (?:) by count. toList() --> wait for all dedup-processes to finish
    flatMap --> Emit every record in the list separatly (for groupTuple)
    groupTuple by readgroup and family
    Map only the best species per family*/
    .map{[it[0],it[1],it[2],it[3],it[4],it[5].text]} //get count-value from text
    .toSortedList({
        a,b -> a[2]+a[1] <=> b[2]+b[1] ?: a[-2] as int <=> b[-2] as int
        })
    .flatMap{n -> n[0..-1]}
    .groupTuple(by:[2,1])   //[[sp,sp,sp], family, rg, [bam,bam,bam],[covered_bp < covered_bp < covered_bp],[count, count, count]]
    .map{n -> [n[0][-1], n[1], n[2], n[3][-1],n[4][-1], n[5][-1]]} //[species, family, rg, bamfile, covered_bp, count]
    .branch{
        no_bed: it[1] in capture_families
        bed: true
    }
    .set{best_deduped}


Channel.fromPath("${params.bedfiles}/*.bed", type:'file')       //all the bedfiles
    .map{[it.baseName.replaceAll(/.masked/,""), it] }           //make a map, with species as key
    .cross(best_deduped.bed)                                    //throw it together
    .map{x, y -> [y[0],y[1],y[2],y[3],x[1],y[4]]}               //get species, family, rg, bamfile and bedfile, covered_bp, (count is ignored here)
    .set{to_bed}

//and filter out reads that intersect with masked regions
process runIntersectBed{
    tag "$rg:$family:$species"
    publishDir 'out', mode: 'copy', saveAs: {out_bam}

    input:
    set species, family, rg, "inbam.bam", "inbed.bed", coverage from to_bed

    output:
    file "outbam.bam"
    set family, rg, species, stdout into (bedfilter_count, bedfilter_data)
    set family, rg, species, "outbam.bam", stdout, coverage into deam_stats
    
    script:
    if(params.byrg){
        out_bam = "${rg}/bed/${family}.${species}_deduped_bedfiltered.bam"
    } else {
        out_bam = "${family}/bed/${rg}.${species}_deduped_bedfiltered.bam"
    }
    """
    bedtools intersect -a inbam.bam -b inbed.bed -v > outbam.bam
    samtools view -c outbam.bam
    """
    }
    
bedfilter_count
    .collectFile(storeDir: 'stats') { family, rg, species, count ->
        [ "${rg}_bedfiltered.tsv", "${family}\t${species}\t${count}"]
    }

//get the captured_families again
best_deduped.no_bed
    .map{[it[1], it[2], it[0], it[3], it[5].trim(), it[4].trim()]}
    .into{captured_no_bed;bedfilter_data_placeholder}

deam_stats
    .concat(captured_no_bed)
    .map{[it[1], it[0], it[2], it[3], it[4].trim().toInteger(), it[5].trim()]}
    .into{ deam_stats;total_rg }

//Sum the counts over the reagroups
total_rg.map{[it[0],it[4]]}
    .groupTuple()     //[RG, [count, count, count, ...]]
    .map{[it[0], it[1].sum()]}
    .set{total_rg}

//throw it together again
deam_stats.combine(total_rg, by:0)
    .into{deam_stats;deam_report}


process reportDamage{
    tag "$rg:$family:$species"
    
    input:
    set rg, family, species, "input.bam", count, coverage, total_rg from deam_stats
    
    output:
    set rg, family, species, stdout into (deam_stats_file, deam_stats_data)
    
    when:
    params.analyze
    
    script:
    """
    bam_deam_stats.py input.bam ${rg} ${family} ${species} ${count} ${coverage} ${total_rg}
    """
}

deam_stats_file
    .collectFile(storeDir: 'stats', keepHeader: true){ rg, family, species, content -> [
        "${rg}_deamination_stats.tsv", content ]
    }


//
//
// Final report
//
//

if(params.report){

//In case the bedfiltering is scipped: replace the value with 'NA'
bedfilter_data_placeholder //[species, family, rg, bamfile, covered_bp, count]
    .map{family, rg, species, bamfile, covered_bp, count -> [family, rg, species, 'NA']}
    .set{bedfilter_data_placeholder}

bedfilter_data
    .concat(bedfilter_data_placeholder)
    .set{bedfilter_data}

//In case --analyze is not set, use a different channel from before 
deam_report.map{[it[0],it[1],it[2],it[4] as String]}.set{deam_report}
post_deam = params.analyze ? deam_stats_data : deam_report

//split and fill with 'NA' if the format doesnt fit to the reportDamage output
post_deam
    .branch {
        empty: it[3].split('\n').size() == 1
        correct: true
    }
    .set{post_deam}

post_deam.empty
    .map{[it[1], it[0], it[2], 'NA','NA','NA','NA','NA','NA']}
    .set{post_deam_empty}

post_deam.correct
    .map{it[3].split('\n')[1].split('\t').flatten()}
    .map{rg, ancient, fam, sp, perc, count, cov, deam5, deam3, cond5, cond3 -> [fam,rg, sp, perc, ancient, deam5, deam3, cond5, cond3]}
    .concat(post_deam_empty)
    .set{post_deam}

//Finally, conbine everything!
mapping_data.combine(deduped_data,by:[1,0,2])
    .combine(coverage_data, by:[1,0,2])
    .combine(bedfilter_data, by:[1,0,2])
    .combine(extraction_data, by:[1,0])
    .combine(post_deam, by:[1,0,2])
    .unique()
    .collectFile(seed:'RG\tFamily\tSpecies\tReadsExtracted\tReadsMapped\tReadsDeduped\tCoveredBP\tReadsBedfiltered\tFamPercentage\tAncientness\tDeam5\tDeam3\tCondDeam5\tCondDeam3', 
        storeDir:'.', newLine:true, sort:true){
        fa,rg,sp,map,dd,cov,bed,bam,ids,ex,p,a,d5,d3,cd5,cd3 -> [
                'final_report.tsv', 
                "${rg}\t${fa}\t${sp}\t${ex.trim()}\t${map.trim()}\t${dd.text.trim()}\t${cov.trim()}\t${bed.trim()}\t${p}\t${a}\t${d5}\t${d3}\t${cd5}\t${cd3}"
        ]
    }
}


