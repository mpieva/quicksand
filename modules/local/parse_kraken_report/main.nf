process PARSE_KRAKEN_REPORT{
    container (workflow.containerEngine ? "merszym/bam_deam:nextflow" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"
    tag "$meta.id"

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path("parsed_report.tsv"), emit: parsed_report
    path 'versions.yml'                       , emit: versions

    script:
    """
    parse_report.py ${report} ${params.krakenuniq_min_kmers} ${params.krakenuniq_min_reads} > parsed_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -f2 -d ' ')
    END_VERSIONS
    """
}