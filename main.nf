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
    set rg, "${dedupfile}" into dedup

    script:
    rg = "${input_bam.baseName}"
    dedupfile = "${rg}-dedup.bam"
    """
    countdups.py -o $dedupfile -s stat.txt -c $cutoff $input_bam
    """
}

process toFasta {
    conda "$baseDir/envs/bam2fasta.yaml"
    // publishDir 'data'

    input:
    set rg, 'input.bam' from dedup

    output:
    set rg, "${fastafile}" into dedup_fasta

    script:
    fastafile = "${rg}-dedup.fa"
    """
    bam2fastx -a -Q -A -o ${fastafile} input.bam
    """
}

process runKraken {
    conda "$baseDir/envs/sediment.yaml"
    // publishDir 'data'

    input:
    set rg, 'input.fa' from dedup_fasta

    output:
    set rg, 'out.translate' into (assigned_to_filter, assigned_to_report)

    script:
    // krakenfile = "${rg}.kraken"
    // kraken_translate = "${rg}.translate"
    """
    kraken -db $params.db --output out.kraken input.fa
    kraken-translate -db $params.db --mpa-format out.kraken >out.translate
    """
}

process reportKraken {
    conda "$baseDir/envs/sediment.yaml"
    publishDir 'kraken'

    input:
    set rg, 'kraken.translate' from assigned_to_report

    output:
    file "${krakenreport}"

    script:
    krakenreport = "${rg}.report"
    """
    kraken-mpa-report --db $params.db kraken.translate >$krakenreport
    """
}
