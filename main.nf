#!/usr/bin/env nextflow

params.bam = ''
params.rg = ''
params.db = '/mnt/ramdisk/refseqReleaseKraken'
params.cutoff = 35

bamfile = file(params.bam)
indices = file(params.rg)
cutoff = params.cutoff

process splitBam {
    conda "$baseDir/envs/sediment.yaml"

    input:
    file 'input.bam' from bamfile
    file 'indices.tsv' from indices

    output:
    file 'split/*.bam' into splitfiles mode flatten

    """
    mkdir split
    splitbam -d split -f indices.tsv --minscore 10 --maxnumber 0 input.bam
    """
}

process removeDups {
    conda "$baseDir/envs/sediment.yaml"
    //publishDir 'data'

    input:
    file input_bam from splitfiles
    
    output:
    set rg, 'dedup.bam' into dedup_to_assign, dedup_to_extract

    script:
    rg = "${input_bam.baseName}"
    """
    countdups.py -o dedup.bam -s stat.txt -c $cutoff $input_bam
    """
}

process toFasta {
    conda "$baseDir/envs/bam2fasta.yaml"
    // publishDir 'data'

    input:
    set rg, 'input.bam' from dedup_to_assign

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
    sleep 5
    """
}
