#!/usr/bin/env nextflow

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

// Inhouse paths
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
params.keeppaired     = ""     // keep paired reads
params.filterunmapped = ""     // filter out unmapped (in case of pre-mapping)
params.analyze        = ""
params.report         = ""
params.splitfile      = "stats/splitcounts.tsv"
params.taxlvl         = "f"    // extract reads on family or order level

if(params.taxlvl !in ['f','o']){
    log.info """
    [quicksand]: ${red}ArgumentError: taxlvl must be one of [o, f] not ${params.taxlvl} ${white}
    """.stripIndent()
    exit 0
}

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
    Use: nextflow run mpieva/quicksand {--rg FILE --bam FILE | --split DIR}
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
        [quicksand]: ${yellow}WARNING: The file 'taxid_map.tsv' is missing in your genomes dir!
        In this run you use the fallback ${baseDir}/assets/taxid_map_example.tsv 
        based on the Refseq release 209. This might lead to WRONG OR MISSING assignments. 
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
params.capture = ''
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
        .map{[it.baseName, it]}
        .set{splitfiles}
} else {
    splitfiles
        .map{[it.baseName, it]}      
        .set{splitfiles}
}

//handle fastq-input
splitfiles
    .branch{
        bam: it[1].getExtension() == "bam"
        fastq: has_ending(it[1], ["fastq","fastq.gz"])
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
post_filter_paired = params.keeppaired ? splitfiles : filter_paired_out

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
    .map{rg,fc,sc -> "${rg}\t${sc.trim()}\t${fc.trim()}"}
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

taxid.splitCsv(sep:'\t')
    .map{[it[0], it[2]]} //tax_id, family, species, (order) -> [tax_id, species]}
    .unique()
    .groupTuple()
    .set{taxid}
    
best_species
    .map{rg, parsed_record -> [rg, parsed_record.splitCsv(sep:'\t',skip:1)]}
    .transpose() // [rg, [Fam, order, taxid, reads, kmers, cov, dup]]
    .map{rg, best_record -> [best_record[2], rg, best_record[0], best_record[1],best_record[3..-1]].flatten()}
    .combine(taxid, by:0)  //[taxid, readgroup, family, order, reads,kmers,cov,dup,[species, species]]
    .into{best_species;kmer_report}

kmer_report
    .map{taxid,rg,fam,order,reads,kmers,cov,dup,species -> [order,fam,rg,reads,kmers,cov,dup]}
    .set{kmer_report}

// This block replaces the default mappings assigned by kraken by those
// specified in the specmap file

def famList = []
if (new File("${params.specmap}").exists()){
    new File("${params.specmap}").eachLine{famList << it.split("\t").flatten()[0]}

    specs = Channel.fromPath("${params.specmap}", type: "file")
        .splitCsv(sep:'\t')
        .map{fam, sp -> [sp.split(",").flatten(),fam]} //[[species,species], Family]
 
    best_species
        .map{taxid,rg,fam,order,reads,kmers,cov,dup,species -> [rg+fam,fam,rg+order,order,species]} //[newKey1, family, newKey2, order [species, species]]
        .branch {
            replace: it[1] in famList
            keep: true
        }
        .set{best_species}

    best_species.replace
        .combine(specs, by:1) //[family, newKey1, newKey2, order, [original_species, ...],[new_species,...])]
        .map{fam,key1,key2,order,species,new_species -> [key1, key2, fam, order, new_species]}
        .set{best_species_replace} //[test3Hominidae, test3Primates, Hominidae, Primates, [Homo_sapiens]]

    best_species.keep
        .map{key1,fam,key2,order,species -> [key1,key2,fam,order,species]}
        .mix(best_species_replace)
        .transpose() // [newKey1, newKey2, Family, Order, Species]
        .into{best_species_post_order; best_species_post_family} 

} else {
    best_species
        .transpose()
        .map{taxid,rg,fam,order,reads,kmers,cov,dup,species -> [rg+fam,rg+order,fam,order,species]}
        .into{best_species_post_order; best_species_post_family}  //[key1, key2, family, order, species]
}

best_species_post_order //Put order-specific newKey2 on first position for the .cross() statement later
    .map{key1,key2,fam,order,sp -> [key2,fam,order,sp]}
    .set{best_species_post_order}


for_extraction
    .cross(kraken_assignments)
    .map{extr,ass -> [extr[0], extr[1], ass[1], ass[1].readLines()]}
    .branch{
        assigned_taxa: it[3].any{it =~ /c__Mammalia.*${params.taxlvl}__./}
        empty: true
    }
    .set{for_extraction}

for_extraction.empty
    .view{"[quicksand]: ${yellow}Info: No Kraken-assignments for Readgroup ${it[0]}${white}"}
    .collectFile(storeDir: 'stats') {it -> [ "${it[0]}_extracted.tsv", "\t"]}

for_extraction.assigned_taxa
    .transpose()
    .filter{it[3] =~ /c__Mammalia.*${params.taxlvl}__./}
    .map{rg,bam,kraken,taxon -> [rg, bam, kraken, (taxon =~ /${params.taxlvl}__([^|]*)/)[0][1]]}
    .unique()
    .ifEmpty{error "----\n${white}[quicksand]:${red} WorkflowError: No families assigned by Kraken at all. Check Input and Database! Exit pipeline${white}"}
    .set{for_extraction}

process gatherByTaxon {
    tag "$rg:$taxon"

    input:
    set rg, 'input.bam', 'kraken.translate', taxon from for_extraction

    output:
    set taxon, rg, 'input.bam', 'ids.txt', stdout into (prepared_for_extraction, count_for_stats, extraction_data)

    script:
    """
    grep "c__Mammalia.*${params.taxlvl}__${taxon}" kraken.translate | cut -f1 | sed "s/\\/[12]//" | tee ids.txt | wc -l
    """
}

count_for_stats
        .collectFile(storeDir: 'stats') {taxon, rg, bamf, idf, count ->
            [ "${rg}_extracted.tsv", "${taxon}\t${count}"]
        }

process extractBam {
    publishDir 'out', mode: 'copy', saveAs: {out_bam}
    tag "$rg:$taxon"

    input:
    set taxon, rg, 'input.bam', 'ids.txt', idcount from prepared_for_extraction

    output:
    set rg, taxon, 'output.bam' into extracted_reads

    when:
    idcount.toInteger() >= params.min_reads 

    script:
    if(params.byrg){
        out_bam = "${rg}/${taxon}/${rg}_extractedReads-${taxon}.bam"
    } else {
        out_bam = "${taxon}/${rg}_extractedReads-${taxon}.bam"
    }
    """
    bamfilter -i ids.txt -l $params.level -o output.bam input.bam
    """
}

extracted_reads // [rg, taxon, extracted_bam] --> newKey1 is rg+fam, newKey2 is rg+order
    .map{rg,taxon,extract_bam -> [rg+taxon, rg, taxon, extract_bam]}
    //extracted reads --> [KeyN, rg, taxon, ExtractedReads_XXX.bam]
    //best_species --> [Newkey1, newkey2, family, order, species]
    .branch{
        family: it[2] =~ ".*idae"
        order: true
    }
    .set{extracted_reads}  

extracted_reads.family //[newKey1, gr, family, extractedReads.bam]
    .cross(best_species_post_family)
    .map{extr, best_fam -> [extr[1],best_fam[3], best_fam[2], best_fam[4], extr[3]]} //[rg, order, family, species, ExtractedReads.bam]
    .set{extracted_reads_families}

extracted_reads.order
    .cross(best_species_post_order)
    .map{extr,best_order -> [extr[1], best_order[2], best_order[1], best_order[3], extr[3]]} //[rg, order, family, species, ExtractedReads.bam]
    .mix(extracted_reads_families) // Mix back together
    .set{extracted_reads}

reference = Channel.fromPath("${params.genome}", type:"dir")
extracted_reads.combine(reference).set{pre_bwa}

process mapBwa {
    publishDir 'out', mode: 'copy', saveAs: {out_bam}, pattern: '*.bam'
    tag "$rg:$family:$species"

    input:
    set rg, order, family, species, "input.bam", genomes from pre_bwa

    output:
    set order, family, rg, species, 'output.bam', taxon into mapped_bam
    set order, family, rg, species, stdout into (mapped_count, mapping_data)

    script:
    if(params.taxlvl == 'o'){
    taxon = order
    } else {
    taxon = family
    }
    if(params.byrg){
        out_bam = "${rg}/${taxon}/aligned/${family}.${species}.bam"
    } else {
        out_bam = "${taxon}/aligned/${rg}.${family}.${species}.bam"
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
        .collectFile(storeDir: 'stats') {order, family, rg, species, count ->
            [ "${rg}_mapped.tsv", "${order}\t${family}\t${species}\t${count}"]
        }


process dedupBam {
    publishDir 'out', mode: 'copy', pattern: "*.bam", saveAs: {out_bam}
    tag "$rg:$family:$species"

    input:
    set order, family, rg, species, 'input.bam', taxon from mapped_bam

    output:
    set order, species, family, rg, "output.bam", stdout, "count.txt", taxon into deduped_bam
    set order, family, rg, species, stdout into (coverage_count, coverage_data)
    set order, family, rg, species, "count.txt" into (deduped_count, deduped_data)

    script:
    if(params.byrg){
        out_bam = "${rg}/${taxon}/aligned/${family}.${species}_deduped.bam"
    } else {
        out_bam = "${taxon}/aligned/${rg}.${family}.${species}_deduped.bam"
    }
    """
    $params.bamrmdup -r -o output.bam input.bam > rmdup.txt
    samtools coverage -H output.bam | cut -f 5
    samtools view -c output.bam > count.txt
    """
}

coverage_count
    .collectFile(storeDir: 'stats') {order, family, rg, species, coverage ->
        [ "${rg}_mapped_coverage.tsv", "${order}\t${family}\t${species}\t${coverage}"]
    }

deduped_count
    .collectFile(storeDir: 'stats') {order, family, rg, species, count_file ->
        [ "${rg}_unique_mapped.tsv", "${order}\t${family}\t${species}\t" + count_file.text]
    }

deduped_bam
    /*sort by readgroup, family and if same (?:) by count. toList() --> wait for all dedup-processes to finish
    flatMap --> Emit every record in the list separatly (for groupTuple)
    groupTuple by readgroup and family
    Map only the best species per family*/
    .map{order,sp,fam,rg,dd_bam,cover,dd_count,taxon -> 
        [order,sp,fam,rg,dd_bam,cover.strip(),dd_count.text.strip(),taxon]}
    .toSortedList({
        a,b -> a[3]+a[2] <=> b[3]+b[2] ?: a[-3] as int <=> b[-3] as int
        })
    .flatMap{n -> n[0..-1]}
    .groupTuple(by:[3,2])   //[[order,order,order][sp,sp,sp], family, rg, [bam,bam,bam],[covered_bp < covered_bp < covered_bp],[count, count, count], [taxon,taxon,taxon]]
    .map{n -> [n[1][-1],n[0][-1], n[2], n[3], n[4][-1],n[5][-1], n[6][-1], n[7][-1]]} //[species, order, family, rg, bamfile, covered_bp, count, taxon]
    .branch{
        no_bed: it[2] in capture_families
        bed: true
    }
    .set{best_deduped}


Channel.fromPath("${params.bedfiles}/*.bed", type:'file') //all the bedfiles
    .map{bedfile -> [bedfile.baseName.replaceAll(/.masked/,""), bedfile]}
    .cross(best_deduped.bed)
    .map{bed, dedup -> [dedup[0],dedup[1],dedup[2],dedup[3],dedup[4],bed[1],dedup[5],dedup[7]]}    
    .set{to_bed} //[species, order, family, rg, bamfile, bedfile, covered_bp, taxon]

//and filter out reads that intersect with masked regions
process runIntersectBed{
    tag "$rg:$family:$species"
    publishDir 'out', mode: 'copy', saveAs: {out_bam}

    input:
    set species, order, family, rg, "inbam.bam", "inbed.bed", coverage, taxon from to_bed

    output:
    file "outbam.bam"
    set order,family, rg, species, stdout into (bedfilter_count, bedfilter_data)
    set order,family, rg, species, "outbam.bam", stdout, coverage into deam_stats
    
    script:
    if(params.byrg){
        out_bam = "${rg}/${taxon}/bed/${family}.${species}_deduped_bedfiltered.bam"
    } else {
        out_bam = "${taxon}/bed/${rg}.${family}.${species}_deduped_bedfiltered.bam"
    }
    """
    bedtools intersect -a inbam.bam -b inbed.bed -v > outbam.bam
    samtools view -c outbam.bam
    """
    }
    
bedfilter_count
    .collectFile(storeDir: 'stats') {order, family, rg, species, count ->
        [ "${rg}_bedfiltered.tsv", "${order}\t${family}\t${species}\t${count}"]
    }

//no bedfiltered families --> prepare for merge with deam_stats
best_deduped.no_bed //[species, order, family, rg, bamfile, covered_bp, count, taxon]
    .map{sp,order,fam,rg,dd_bam,cover,dd_count,taxon -> 
        [order,fam,rg,sp,dd_bam,dd_count,cover]}
    .into{captured_no_bed; bedfilter_data_placeholder}

deam_stats //[order, family, rg, species, bamfile, count, coverage]
    .concat(captured_no_bed)
    .map{order,fam,rg,sp,bed_bam,bed_count,cover -> 
        [rg,order,fam,sp,bed_bam,bed_count.toInteger(),cover]}
    .into{deam_stats;total_rg}

//In order to get the percentage, we need the count per rg. Thus, sum the counts over the reagroups
total_rg.map{rg,order,fam,sp,bed_bam,bed_count,cover -> [rg,bed_count]}
    .groupTuple()     //[RG, [count, count, count, ...]]
    .map{rg,bed_count -> [rg, bed_count.sum()]}
    .set{total_rg}

//throw it together again
deam_stats.combine(total_rg, by:0)
    .map{rg,order,family,species,bam,count,coverage,total_rg -> [rg,order,family,species,bam,count,coverage,total_rg==0?0:(count*100/total_rg).trunc(2)]}
    .into{deam_stats;deam_report} //[rg, order, family, species, bamfile, count, coverage, perc]


process reportDamage{
    tag "$rg:$family:$species"
    
    input:
    set rg, order, family, species, "input.bam", count, coverage, perc from deam_stats
    
    output:
    set rg, order, family, species, stdout, perc into (deam_stats_file, deam_stats_data)
    
    when:
    params.analyze
    
    script:
    """
    bam_deam_stats.py input.bam ${rg} ${order} ${family} ${species} ${count} ${coverage} ${perc}
    """
}

deam_stats_file
    .collectFile(storeDir: 'stats', keepHeader: true){rg, order, family, species, content, perc -> [
        "${rg}_deamination_stats.tsv", content ]
    }


//
//
// Final report
//
//

if(params.report){

//In case the bedfiltering is scipped: replace the value with 'NA'
bedfilter_data_placeholder //[order, family, rg, species, bamfile, count, coverage]
    .map{order, family, rg, species, bamfile, count, coverage -> [order, family, rg, species, 'NA']}
    .set{bedfilter_data_placeholder}

bedfilter_data //values of the actual bedfiltered bamfiles --> [order,family, rg, species, count]
    .concat(bedfilter_data_placeholder)
    .set{bedfilter_data}

//In case --analyze is not set, use a different channel from before damage-analysis
//deam report --> [rg, order, family, species, bamfile, count, coverage, total_count, perc]
deam_report
    .map{rg,order,fam,sp,bed_bam,bed_count,cover,perc -> 
        [rg,order,fam,sp,bed_count as String,perc]}
    .set{deam_report}
post_deam = params.analyze ? deam_stats_data : deam_report

//split and fill with 'NA' if the format doesnt fit to the reportDamage output
post_deam //[rg, order, family, species, count or deam_output, perc]
    .branch {
        empty: it[4].split('\n').size() == 1
        correct: true
    }
    .set{post_deam}

post_deam.empty //[rg, order, family, species, count, perc]
    .map{rg,order,fam,sp,bed_count,perc -> [order,fam,rg,sp,perc,'NA','NA','NA','NA','NA']}
    .set{post_deam_empty} //[order,fam,rg,sp,perc,ancient,deam5,deam3,cond5,cond3]

post_deam.correct
    .map{rg,order,fam,sp,deam_out,perc -> [deam_out.split('\n')[1].split('\t')].flatten()}
    .map{rg,ancient,order,fam,sp,perc,count,cov,deam5,deam3,cond5,cond3 -> 
        [order,fam,rg,sp,perc,ancient,deam5,deam3,cond5,cond3]}
    .concat(post_deam_empty) //mix together with non-analyzed
    .set{post_deam}

//Finally, conbine everything!
//Mapping_data: [order, family, rg, species, mapping_count]
mapping_data
    .combine(deduped_data,by:[0,2,1,3])
    .combine(coverage_data, by:[0,2,1,3])
    .combine(bedfilter_data, by:[0,2,1,3])
    .set{mapping_data}

//To combine the extraction data, we have to use either the order or the family...
if(params.taxlvl == 'f'){
    mapping_data
        .map{order,fam,rg,sp,map,ded,cov,bed -> [fam,rg,order,sp,map,ded,cov,bed]}
        .combine(extraction_data, by:[0,1])
	.map{fam,rg,order,sp,map,ded,cov,bed,bam,ids,ex -> [fam,rg,order,sp,map.trim() as int,ded.text.trim() as int,cov.trim() as int,bed.trim() as int,bam,ids,ex.trim() as int]}
 	.map{fam,rg,order,sp,map,ded,cov,bed,bam,ids,ex -> [order,fam,rg,sp,map,ded,cov,bed,bam,ids,ex,(ex==0 || map==0) ? 0 : map/ex, (ded == 0 || map==0) ? 0: map/ded]}
        .set{mapping_data}
}else{
    mapping_data
        .map{order, fam, rg, sp, map, ded, cov, bed -> [order,rg,fam,sp,map,ded,cov,bed]}
        .combine(extraction_data, by:[0,1])
	.map{order,rg,fam,sp,map,ded,cov,bed,bam,ids,ex -> [fam,rg,order,sp,map.trim() as int,ded.text.trim() as int,cov.trim() as int,bed.trim() as int,bam,ids,ex.trim() as int]}
        .map{fam,rg,order,sp,map,ded,cov,bed,bam,ids,ex -> [order,fam,rg,sp,map,ded,cov,bed,bam,ids,ex,(ex==0 || map==0) ? 0 : map/ex, (ded == 0 || map==0) ? 0: map/ded]}
        .set{mapping_data}
}
mapping_data
    .combine(post_deam, by:[0,2,1,3])
    .combine(kmer_report, by:[0,2,1])
    .collectFile(seed:'RG\tFamilyReads\tFamilyKmers\tKmerCoverage\tKmerDupRate\tOrder\tFamily\tSpecies\tExtractLVL\tReadsExtracted\tReadsMapped\tProportionMapped\tReadsDeduped\tDuplicationRate\tCoveredBP\tReadsBedfiltered\tFamPercentage\tAncientness\tDeam5\tDeam3\tCondDeam5\tCondDeam3', 
        storeDir:'.', newLine:true, sort:true){
        order,fa,rg,sp,map,dd,cov,bed,bam,ids,ex,mapdrop,dup,p,a,d5,d3,cd5,cd3,read,kmer,kmercov,kmerdup -> [
                'final_report.tsv', 
                "${rg}\t${read}\t${kmer}\t${kmercov}\t${kmerdup}\t${order}\t${fa}\t${sp}\t${params.taxlvl}\t${ex}\t${map}\t${mapdrop}\t${dd}\t${dup}\t${cov}\t${bed}\t${p}\t${a}\t${d5}\t${d3}\t${cd5}\t${cd3}"
        ]
    }
}


