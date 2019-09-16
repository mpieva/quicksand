#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    Usage: nextflow run sediment_nf --bam PATH --rg PATH --db PATH --genome PATH
                                    [--cutoff N] [--level N] [--filter] [--dedup]
           
    Run sediment analysis pipeline.
    
    required arguments:
      --bam PATH      input BAM file
      --rg  PATH      tab-separated file containing index combinations
                      - format of file: 'LibID<tab>P7<tab>P5'
      --db PATH       Kraken database
      --genome PATH   genome for alignment
      
    optional arguments:
      --cutoff N      length cutoff (default: 35)
      --filter        filter out unmapped and paired reads
      --dedup         dedupicate input
      --level N       set BGZF compression level (default: 6)
      
    A selection of built-in Nextflow flags that may be of use:
      -resume         resume processing; do not re-run completed processes
      -profile NAME   use named execution profile
      -qs N           queue size; number of CPU cores to use (default: all)
      -N EMAIL        send completion notifcation to email address
    
    The following execution profiles are available (see '-profile' above):
      * standard      execute all processes on local host
      * cluster       execute certain CPU-intenisve processes on SGE cluster
    """.stripIndent()
}

params.help = false
if (params.help) {
    helpMessage()
    exit 0
}

params.bam    = ''
params.rg     = ''
params.db     = ''
params.genome = ''
params.cutoff = 35
params.filter = false       // filter out unmapped and paired
params.dedup  = false       // deduplicate after splitting
params.level  = 0           // bgzf compression level for intermediate files, 0..9
params.bwa    = '/home/public/usr/bin/bwa'  // not intended to be set by user


process splitBam {
    conda "$baseDir/envs/sediment.yaml"
    maxForks 1
    label 'local'

    input:
    file 'input.bam' from file(params.bam)
    file 'indices.tsv' from file(params.rg)

    output:
    file '*.bam' into splitfiles mode flatten

    script:
    """
    splitbam -c $params.level -f indices.tsv --minscore 10 --maxnumber 0 input.bam
    """
}

splitfiles
    .map { [it.baseName, it] }
    .set { splitfiles }

filter_in = params.filter ? splitfiles : Channel.empty()

process filterBam {
    conda "$baseDir/envs/sediment.yaml"
    tag "$rg"

    input:
    set rg, 'input.bam' from filter_in

    output:
    set rg, 'output.bam' into filter_out

    script:
    """
    samtools view -b -u -F 5 -o output.bam input.bam
    """
}

post_filter = params.filter ? filter_out : splitfiles

dedup_in = params.dedup ? post_filter : Channel.empty()

process removeDups {
    conda "$baseDir/envs/sediment.yaml"
    tag "$rg"

    input:
    set rg, 'input.bam' from dedup_in
    
    output:
    set rg, 'dedup.bam' into dedup_out

    script:
    """
    countdups.py -l $params.level -o dedup.bam -s stat.txt -c $params.cutoff input.bam
    """
}

(params.dedup ? dedup_out : post_filter)
    .into{ tofasta_in; for_extraction }

process toFasta {
    conda "$baseDir/envs/bam2fasta.yaml"
    tag "$rg"

    input:
    set rg, 'input.bam' from tofasta_in

    output:
    set rg, 'output.fa' into dedup_fasta

    script:
    """
    bam2fastx -a -Q -A -o output.fa input.bam
    """
}

process runKraken {
    conda "$baseDir/envs/sediment.yaml"
    publishDir 'kraken', mode: 'link'
    memory '16GB'   // XXX make dynamic based on size of params.db
    label 'bigmem'
    label 'local'   // Currently having issued with bigmem jobs on SGE
    tag "$rg"

    input:
    set rg, 'input.fa' from dedup_fasta

    output:
    set rg, "${kraken_translate}" into kraken_assignments
    set rg, "${kraken_out}" into kraken_raw

    script:
    kraken_out = "${rg}.kraken"
    kraken_translate = "${rg}.translate"
    """
    kraken --db $params.db --output $kraken_out input.fa
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

//for_extraction = ( params.dedup ? dedup_to_extract : split_to_extract )

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

process extractBam {
    conda "$baseDir/envs/sediment.yaml"
    publishDir 'out/blast', mode: 'link', saveAs: { out_bam }
    tag "$rg:$family"

    input:
    set rg, 'input.bam', 'kraken.translate', family from for_extraction

    output:
    set family, rg, 'output.bam' into extracted_reads
    set family, rg, stdout into extracted_read_count

    script:
    out_bam = "${family}/${rg}_extractedReads-${family}.bam"
    """
    extract_bam.py -f $family -k kraken.translate -c $params.level -o output.bam input.bam
    """
}

extracted_read_count
    .collectFile(storeDir: 'out/blast') { family, rg, count ->
        [ "${rg}.tsv", "${family}\t${count}"]
    }

Channel.fromPath("${params.genome}/*", type: 'dir')
    .map { [it.baseName, file("${it}/*.fasta")] }
    .cross(extracted_reads)
    .map { x, y -> [x[0], y[1], y[2], x[1]] }   // family, rg, bamfile, fasta_list
    .transpose(by: 3)
    .set { for_mapping }

process mapBwa {
    conda "$baseDir/envs/sediment.yaml"
    publishDir 'out/blast', mode: 'link', saveAs: { out_bam }
    tag "$rg:$family"

    input:
    set family, rg, 'input.bam', genome_fasta from for_mapping

    output:
    file 'output.bam'

    script:
    species = genome_fasta.baseName
    out_bam = "${family}/aligned/${rg}.${species}.bam"
    """
    $params.bwa bam2bam -g $genome_fasta -n 0.01 -o 2 -l 16500 --only-aligned input.bam \
    | samtools sort -o output.bam
    """
}
