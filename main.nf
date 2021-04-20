#!/usr/bin/env nextflow

red = "\033[0;31m"
white = "\033[0m"
cyan = "\033[0;36m"
yellow = "\033[0;33m"

log.info """
[Sediment_nf]: Execution started: ${workflow.start.format('dd.MM.yyyy HH:mm')} ${cyan}

  _____         _  _                   _      _____  _____ 
 |   __| ___  _| ||_| _____  ___  ___ | |_   |   | ||   __|
 |__   || -_|| . || ||     || -_||   ||  _|  | | | ||   __|
 |_____||___||___||_||_|_|_||___||_|_||_|    |_|___||__|   
 ${white}${workflow.manifest.description} ${cyan}~ Version ${workflow.manifest.version} ${white}

==============================================================
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

// User-independent params
params.cutoff         = 35
params.quality        = 25
params.bwacutoff      = 0
params.krakenthreads  = 4      // number of threads per kraken process
params.level          = 0      // bgzf compression level for intermediate files, 0..9
params.keeppaired     = false  // keep paired reads
params.filterunmapped = false  // filter out unmapped
params.krakenfilter   = false  // a kraken-filter step with the given weight
params.analyze        = ""
params.splitfile      = "stats/splitcounts.tsv"

def env = System.getenv()

def validate_dir(path, flag){
    File test = new File("$path")
    if (!(test.exists() && test.isDirectory())){
        log.info "[sediment_nf]: ${red}InputError: The $flag directory doesn't exist ${white}\nPath: $path"
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
    [sediment_nf]: ${red}ArgumentError: Too many arguments ${white}
    Use: nextflow run sediment_nf {--rg FILE --bam FILE | --split DIR }
    """.stripIndent()
    exit 0
} 
if(!params.split && !(params.bam && params.rg)){
    log.info """
    [sediment_nf]: ${red}ArgumentError: Too few arguments ${white}
    See: nextflow run sediment_nf --help
    """.stripIndent()
    exit 0
}

inbam = params.bam? Channel.fromPath("${params.bam}") : Channel.empty()
indexfile = params.rg? Channel.fromPath("${params.rg}") : Channel.empty()


// Validate: GENOME and TAXIDMAP
params.genome = env["SED_GENOME"]
if(!params.genome){
    log.info """[sediment_nf]: ${red}ArgumentError: Missing --genome flag ${white}"""
    exit 0;
} else { 
    validate_dir("${params.genome}","--genome") 
}

if(new File("${params.genome}/taxid_map.tsv").exists()){
    taxid = Channel.fromPath("${params.genome}/taxid_map.tsv", type:'file')
} else {
    taxid = Channel.fromPath("${baseDir}/assets/taxid_map_example.tsv", type:'file')
    log.info """
        [sediment_nf]: ${yellow}CRITICAL WARNING: The file 'taxid_map.tsv' is missing in your genomes dir!
        In this run you use the fallback ${baseDir}/assets/taxid_map_example.tsv 
        based on the Refseq release 202. This might lead to WRONG OR MISSING assignments. 
        Please check the Docs, how to set up an own taxid_map.tsv ${white}
        """.stripIndent()
}

// Validate: KRAKEN DB
params.db = env["SED_DB"]
if(!params.db){
    log.info """[sediment_nf]: ${red}ArgumentError: Missing --db flag ${white}"""
    exit 0;
} else {
    validate_dir("${params.db}", "--db")
}


// Validate: BEDFILES
params.bedfiles = env["SED_BEDFILES"]
if(!params.bedfiles){
    log.info """[sediment_nf]: ${red}ArgumentError: Missing --bedfiles flag ${white}"""
    exit 0;
} else {
    validate_dir("${params.bedfiles}", "--bedfiles")
}

// Optional
params.specmap = env["SED_SPECMAP"]

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
        log.info "[sediment_nf]: ${red}InputError: The --split directory doesn't exist ${white}\nPath: ${params.split}"
        exit 0;
    }
    Channel.fromPath("${params.split}/*")
        .ifEmpty{error "[sediment_nf]: InputError: The Split-Directory is empty"}
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
    .view{"[sediment_nf]: ${yellow}WARNING: ${it[1]} omitted. File has neither bam nor fastq-ending!${white}"}

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
    .ifEmpty{error "----\n${white}[sediment_nf]:${red}WorkflowError: No input-files. Scheck SPLIT-dir or RG-BAM combination. Exit pipeline${white}"}
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

process runKraken {
    cpus "${params.krakenthreads}"
    memory '12GB'
    label 'bigmem'
    label 'local'
    tag "$rg"

    input:
    set rg, "input.fa", db from pre_kraken

    output:
    set rg, "output.kraken", db into kraken_out

    script:
    """
    kraken --threads ${task.cpus} --db ${db} --output output.kraken --fasta-input input.fa
    """
}

kraken_filter_in = params.krakenfilter ? kraken_out : Channel.empty()

process filterKraken{
    publishDir 'kraken', mode: 'copy', saveAs: {"${rg}.kraken_filter"}
    label 'local'
    tag "$rg"

    input:
    set rg, "input.kraken", db from kraken_filter_in

    output:
    set rg, "output_filtered.kraken", db into kraken_filter_out

    when:
    params.krakenfilter

    script:
    """
    kraken-filter -threshold $params.krakenfilter --db ${db} input.kraken > output_filtered.kraken
    """
}

post_kraken_filter = params.krakenfilter ? kraken_filter_out : kraken_out

process translateKraken{
    publishDir 'kraken', mode: 'copy', pattern:"*translate", saveAs: {"${rg}.translate"}
    publishDir 'kraken', mode: 'copy', pattern:"*report", saveAs: {"${rg}.report"}
    label 'local'
    tag "$rg"

    input:
    set rg, "input.kraken", db from post_kraken_filter

    output:
    set rg, "kraken.translate" into kraken_assignments
    set rg, "kraken.report" into find_best

    script:
    """
    kraken-translate --db ${db} --mpa-format input.kraken > kraken.translate
    kraken-report --db ${db} input.kraken > kraken.report
    """
}
process findBestSpecies{
    tag "$rg"

    input:
    set rg, "kraken.report" from find_best

    output:
    set rg, "parsed_record.tsv" into best_species

    script:
    """
    parse_report.py kraken.report
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
    .view{"[sediment_nf]: ${yellow}Info: No Kraken-assignments for Readgroup ${it[0]}${white}"}
    .collectFile(storeDir: 'stats') { it -> [ "${it[0]}_extracted.tsv", "\t"] }

for_extraction.assigned_taxa
    .transpose()
    .filter { it[3] =~ /c__Mammalia.*f__./ }
    .map { rg, bam, kraken, asn -> [rg, bam, kraken, (asn =~ /f__([^|]*)/)[0][1]] }
    .unique()
    .ifEmpty{error "----\n${white}[sediment_nf]:${red} WorkflowError: No families assigned by Kraken at all. Check Input and Database! Exit pipeline${white}"}
    .set { for_extraction }

process gatherByFamily {
    tag "$rg:$family"

    input:
    set rg, 'input.bam', 'kraken.translate', family from for_extraction

    output:
    set family, rg, 'input.bam', 'ids.txt', stdout into (prepared_for_extraction, count_for_stats)

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
    publishDir 'out', mode: 'copy', saveAs: {"${family}/${rg}_extractedReads-${family}.bam"}
    tag "$rg:$family"

    input:
    set family, rg, 'input.bam', 'ids.txt', idcount from prepared_for_extraction

    output:
    set rg, family, 'output.bam' into extracted_reads

    when:
    idcount.toInteger() >= params.bwacutoff

    script:
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
    set family, rg, species, stdout into mapped_count

    script:
    out_bam = "${family}/aligned/${rg}.${species}.bam"
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
    publishDir 'out', mode: 'copy', pattern: "*.bam", saveAs: {"${family}/aligned/${rg}.${species}_deduped.bam"}
    tag "$rg:$family:$species"

    input:
    set family, rg, species, 'input.bam' from mapped_bam

    output:
    set species, family, rg, "output.bam", stdout into deduped_bam
    set family, rg, species, stdout into coverage_count
    set family, rg, species, "count.txt" into deduped_count

    script:
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
    .toSortedList({
        a,b -> a[2]+a[1] <=> b[2]+b[1] ?: a[-1] as int <=> b[-1] as int
        })
    .flatMap{n -> n[0..-1]}
    .groupTuple(by:[2,1])   //[[sp,sp,sp], family, rg, [bam,bam,bam],[count < count < count]]
    .map{n -> [n[0][-1], n[1], n[2], n[3][-1],n[4][-1]]} //[species, family, rg, bamfile, count, covered_bp]
    .set{best_deduped}


Channel.fromPath("${params.bedfiles}/*.bed", type:'file')   //all the bedfiles
    .map{[it.baseName.replaceAll(/.masked/,""), it] }       //make a map, with species as key
    .cross(best_deduped)                                    //throw it together
    .map{x, y -> [y[0],y[1],y[2],y[3],x[1],y[4]]}           //get species, family, rg, bamfile and bedfile, covered_bp
    .set{to_bed}

//and filter out reads that intersect with masked regions
process runIntersectBed{
    tag "$rg:$family:$species"
    publishDir 'out', mode: 'copy', saveAs: {"${family}/bed/${rg}.${species}_deduped_bedfiltered.bam"}

    input:
    set species, family, rg, "inbam.bam", "inbed.bed", coverage from to_bed

    output:
    file "outbam.bam"
    set family, rg, species, stdout into bedfilter_count
    set family, rg, species, "outbam.bam", stdout, coverage into deam_stats
    
    script:
    """
    bedtools intersect -a inbam.bam -b inbed.bed -v > outbam.bam
    samtools view -c outbam.bam
    """
    }
    
bedfilter_count
    .collectFile(storeDir: 'stats') { family, rg, species, count ->
        [ "${rg}_bedfiltered.tsv", "${family}\t${species}\t${count}"]
    }

deam_stats.map{[it[1], it[0], it[2], it[3], it[4].trim().toInteger(), it[5].trim()]}
    .into{ deam_stats;total_rg }

//Sum the counts over the reagroups
total_rg.map{[it[0],it[4]]}
    .groupTuple()     //[RG, [count, count, count, ...]]
    .map{[it[0], it[1].sum()]}
    .set{total_rg}

//throw it together again
deam_stats.combine(total_rg, by:0)
    .map{x -> [x[0],x[1],x[2],x[3],x[4],x[5],x[6]]}
    .set{deam_stats}


process reportDamage{
    tag "$rg:$family:$species"
    
    input:
    set rg, family, species, "input.bam", count, coverage, total_rg from deam_stats
    
    output:
    set rg, stdout into deam_stats_file
    
    when:
    params.analyze
    
    script:
    """
    bam_deam_stats.py input.bam ${rg} ${family} ${species} ${count} ${coverage} ${total_rg}
    """
}

deam_stats_file
    .collectFile(storeDir: 'stats', keepHeader: true){ rg, content -> [
        "${rg}_deamination_stats.tsv", content ]
    }

