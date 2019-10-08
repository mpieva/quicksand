#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    Usage: nextflow run sediment_nf --bam PATH --rg PATH --db PATH --genome PATH
                                    [--cutoff N] [--level N] [--filter] [--dedup]
           
    Run sediment analysis pipeline.
    
    required arguments:
      --bam PATH         input BAM file
      --rg  PATH         tab-separated file containing index combinations
                         - format of file: 'LibID<tab>P7<tab>P5'
      --db PATH          Kraken database
      --genome PATH      genome for alignment
      
    optional arguments:
      --cutoff N         length cutoff (default: 35)
      --quality N        quality filter (default: 25)
      --bwacutoff N      cutoff for number of kraken matches (default: 0)
      --filterpaired     filter paired reads
      --filterunmapped   filter unmapped reads
      --level N          set BGZF compression level (default: 6)
      --krakenthreads N  numbrer of threads per Kraken process (default: 4)
      
    A selection of built-in Nextflow flags that may be of use:
      -resume            resume processing; do not re-run completed processes
      -profile NAME      use named execution profile
      -qs N              queue size; number of CPU cores to use (default: all)
      -N EMAIL           send completion notifcation to email address
    
    The following execution profiles are available (see '-profile' above):
      * standard         execute all processes on local host
      * cluster          execute certain CPU-intenisve processes on SGE cluster
    """.stripIndent()
}

params.help = false
if (params.help) {
    helpMessage()
    exit 0
}

params.bam            = ''
params.rg             = ''
params.db             = ''
params.genome         = ''
params.cutoff         = 35
params.quality        = 25
params.bwacutoff      = 0
params.filterpaired   = false  // filter out paired
params.filterunmapped = false  // filter out unmapped
params.dedup          = false  // deduplicate after splitting
params.level          = 0      // bgzf compression level for intermediate files, 0..9
params.krakenthreads  = 4      // number of threads per kraken process

// The following parameters are not meant to be set by the end user:
params.bwa            = '/home/public/usr/bin/bwa'
params.bamrmdup       = '/home/bioinf/usr/bin/bam-rmdup'


process splitBam {
    conda "$baseDir/envs/sediment.yaml"
    maxForks 1
    publishDir 'split', mode: 'link'
    label 'local'

    input:
    file 'input.bam' from file(params.bam)
    file 'indices.tsv' from file(params.rg)

    output:
    file '*.bam' into splitfiles mode flatten
    stdout into splitscriptstats

    script:
    """
    splitbam -s -c $params.level -f indices.tsv --minscore 10 --maxnumber 0 input.bam
    """
}

splitscriptstats
    .collectFile(storeDir: 'stats', name: "splitstats.tsv", newLine: true)

splitfiles
    .map { [it.baseName, it] }
    .into { splitfiles; splitstats }

process splitStats {
    conda "$baseDir/envs/sediment.yaml"
    tag "$rg"

    input:
    set rg, 'input.bam' from splitstats

    output:
    set rg, stdout into splitcounts

    script:
    """
    samtools view -c input.bam
    """
}

filter_paired_in = params.filterpaired ? splitfiles : Channel.empty()

process filterPaired {
    conda "$baseDir/envs/sediment.yaml"
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

post_filter_paired = params.filterpaired ? filter_paired_out : splitfiles

filter_unampped_in = params.filterunmapped ? post_filter_paired : Channel.empty()

process filterUnmapped {
    conda "$baseDir/envs/sediment.yaml"
    tag "$rg"

    input:
    set rg, 'input.bam' from filter_unampped_in

    output:
    set rg, 'output.bam' into filter_unmapped_out

    script:
    """
    samtools view -b -u -F 4 -o output.bam input.bam
    """
}

post_filter_unmapped = params.filterunmapped ? filter_unmapped_out : post_filter_paired

process filterLength {
    conda "$baseDir/envs/sediment.yaml"
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

filtercounts.join(splitcounts)
    .map { rg, fc, sc -> "${rg}\t${sc.trim()}\t${fc.trim()}"}
    .collectFile(storeDir: 'stats', name: "splitcounts.tsv", newLine: true,
                 seed: "readgroup\tsplit count\tfiltered count")

process toFasta {
    conda "$baseDir/envs/sediment.yaml"
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

process runKraken {
    conda "$baseDir/envs/sediment.yaml"
    publishDir 'kraken', mode: 'link'
    cpus "${params.krakenthreads}"
    memory '16GB'   // XXX make dynamic based on size of params.db
    label 'bigmem'
    label 'local'   // Currently having issues with bigmem jobs on SGE
    tag "$rg"

    input:
    set rg, 'input.fa' from tofasta_out

    output:
    set rg, "${kraken_translate}" into kraken_assignments
    set rg, "${kraken_out}" into kraken_raw

    script:
    kraken_out = "${rg}.kraken"
    kraken_translate = "${rg}.translate"
    """
    kraken --threads ${task.cpus} --db $params.db --output $kraken_out --fasta-input input.fa
    kraken-translate --db $params.db --mpa-format $kraken_out >$kraken_translate
    """
}

process statsKraken {
    conda "$baseDir/envs/sediment.yaml"
    publishDir 'kraken', mode: 'link', saveAs: { "${rg}.report" }
    label 'local'
    tag "$rg"

    input:
    set rg, 'input.kraken' from kraken_raw

    output:
    file 'kraken.report'

    script:
    """
    kraken-mpa-report --db $params.db input.kraken >kraken.report
    """
}

// - Match split/deduped bam to kraken output based on readgroup
// - Filter for families under Mammalia
// - Emit unique (rg, bamfile, kraken_translate_file, family_name)
//
for_extraction
    .cross(kraken_assignments)
    .map { [it[0][0], it[0][1], it[1][1], it[1][1].readLines()] }
    .transpose()
    .filter { it[3] =~ /c__Mammalia.*f__/ }
    .map { rg, bam, kraken, asn -> [rg, bam, kraken, (asn =~ /f__([^|]*)/)[0][1]] }
    .unique()
    .set { for_extraction }

process gatherByFamily {
    tag "$rg:$family"

    input:
    set rg, 'input.bam', 'kraken.translate', family from for_extraction

    output:
    set family, rg, 'input.bam', 'ids.txt', stdout into prepared_for_extraction

    script:
    out_bam = "${family}/${rg}_extractedReads-${family}.bam"
    """
    grep "c__Mammalia.*f__$family" kraken.translate | cut -f1 | tee ids.txt | wc -l
    """
}

prepared_for_extraction
    .into { count_for_stats; prepared_for_extraction }

count_for_stats
        .collectFile(storeDir: 'stats') { family, rg, bamf, idf, count ->
            [ "${rg}_extracted.tsv", "${family}\t${count}"]
        }

process extractBam {
    conda "$baseDir/envs/sediment.yaml"
    publishDir 'out', mode: 'link', saveAs: { out_bam }
    tag "$rg:$family"

    input:
    set family, rg, 'input.bam', 'ids.txt', idcount from prepared_for_extraction

    output:
    set family, rg, 'output.bam' into extracted_reads

    when:
    idcount.toInteger() >= params.bwacutoff

    script:
    out_bam = "${family}/${rg}_extractedReads-${family}.bam"
    """
    bamfilter -i ids.txt -l $params.level -o output.bam input.bam
    """
}

Channel.fromPath("${params.genome}/*", type: 'dir')
    .map { [it.baseName, file("${it}/*.fasta")] }
    .cross(extracted_reads)
    .map { x, y -> [x[0], y[1], y[2], x[1]] }   // family, rg, bamfile, fasta_list
    .transpose(by: 3)
    .set { for_mapping }

process mapBwa {
    conda "$baseDir/envs/sediment.yaml"
    tag "$rg:$family:$species"

    input:
    set family, rg, 'input.bam', genome_fasta from for_mapping

    output:
    set family, rg, species, 'output.bam' into mapped_bam
    set family, rg, species, stdout into mapped_count

    script:
    species = genome_fasta.baseName
    out_bam = "${family}/aligned/${rg}.${species}.bam"
    """
    samtools sort -n -l0 input.bam \
    | $params.bwa bam2bam -g $genome_fasta -n 0.01 -o 2 -l 16500 --only-aligned - \
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
    publishDir 'out', mode: 'link', saveAs: { out_bam }
    tag "$rg:family:$species"

    input:
    set family, rg, species, 'input.bam' from mapped_bam

    output:
    set family, rg, species, stdout into deduped_count
    file 'output.bam'

    script:
    out_bam = "${family}/aligned/${rg}.${species}.bam"
    """
    $params.bamrmdup -r -o output.bam input.bam >rmdup.txt
    samtools view -c output.bam
    """
}

deduped_count
        .collectFile(storeDir: 'stats') { family, rg, species, count ->
            [ "${rg}_unique_mapped.tsv", "${family}\t${species}\t${count}"]
        }
