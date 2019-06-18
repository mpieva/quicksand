#!/usr/bin/env nextflow

params.bam = ''
params.rg = ''
params.db = '/mnt/ramdisk/refseqReleaseKraken'
params.cutoff = 35
params.dedup = false
params.level = 0

bamfile = file(params.bam)
indices = file(params.rg)
cutoff = params.cutoff
level = params.level        // bgzf compression level for intermediate files, 0..9

process splitBam {
    conda "$baseDir/envs/sediment.yaml"

    input:
    file 'input.bam' from bamfile
    file 'indices.tsv' from indices

    output:
    file 'split/*.bam' into splitfiles, split_to_extract mode flatten

    """
    mkdir split
    splitbam -c $level -d split -f indices.tsv --minscore 10 --maxnumber 0 input.bam
    """
}

splitfiles
    .map { [it.baseName, it] }
    .set { splitfiles }

split_to_extract
    .map { [it.baseName, it] }
    .set { split_to_extract }

(dedup_input, tofasta_input) = ( params.dedup
                             ? [splitfiles, Channel.empty()]
                             : [Channel.empty(), splitfiles] )

process removeDups {
    conda "$baseDir/envs/sediment.yaml"
    //publishDir 'data'

    input:
    set rg, 'input.bam' from dedup_input
    
    output:
    set rg, 'dedup.bam' into dedup_to_assign, dedup_to_extract

    script:
    // rg = "${input_bam.baseName}"
    """
    countdups.py -l $level -o dedup.bam -s stat.txt -c $cutoff input.bam
    """
}

process toFasta {
    conda "$baseDir/envs/bam2fasta.yaml"
    // publishDir 'data'

    input:
    set rg, 'input.bam' from tofasta_input.mix(dedup_to_assign)

    output:
    set rg, 'output.fa' into dedup_fasta

    script:
    """
    bam2fastx -a -Q -A -o output.fa input.bam
    """
}

process runKraken {
    conda "$baseDir/envs/sediment.yaml"
    publishDir 'kraken', mode: 'copy'

    input:
    set rg, 'input.fa' from dedup_fasta

    output:
    set rg, "${kraken_translate}" into kraken_assignments
    file "${kraken_out}" into kraken_raw
    file "${kraken_report}" into kraken_stats

    script:
    kraken_out = "${rg}.kraken"
    kraken_translate = "${rg}.translate"
    kraken_report = "${rg}.report"
    """
    kraken -db $params.db --output $kraken_out input.fa
    kraken-translate -db $params.db --mpa-format $kraken_out >$kraken_translate
    kraken-mpa-report --db $params.db $kraken_out >$kraken_report
    """
}
// look into executor.exitReadTimeOut instead of `sleep`

for_extraction = ( params.dedup ? dedup_to_extract : split_to_extract )

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
    publishDir 'out', mode: 'copy', overwrite: true, saveAs: { "${out_bam}" }

    input:
    set rg, 'input.bam', 'kraken.translate', family from for_extraction

    output:
    file 'output.bam' into extracted_reads

    script:
    out_bam = "${rg}_extracted_reads-${family}.bam"
    """
    extract_bam.py -f $family -k kraken.translate -c $level -o output.bam input.bam
    """
}
