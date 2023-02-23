process SAMTOOLS_FILTER {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    label 'local'
    tag "$meta.id:Flag:${params.bamfilterflag}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("filtered_*.bam")  , emit: bam
    tuple val(meta), stdout                  , emit: prefilter_count
    tuple val(meta), path("filtercount.txt") , emit: txt


    script:
    """
    samtools view -c -F 128 ${bam}
    samtools view -b -u -F ${params.bamfilterflag} -o filtered_${bam} ${bam}
    samtools view -c filtered_${bam} > filtercount.txt
    """
}