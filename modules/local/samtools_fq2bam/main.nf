process SAMTOOLS_FQ2BAM{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    tag "$meta.RG"
    label 'local'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.bam"), emit: bam

    script:
    def name = fastq.baseName
    """
    samtools import -0 ${fastq} -o \"${name}.bam\"
    """
}