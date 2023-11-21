process SAMTOOLS_MPILEUP{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    tag "${meta.id}:${meta.Taxon}:${meta.Species}"
    label 'local'

    input:
    tuple val(meta), path(bam1), path(bam3), path(bam)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv

    script:
    def args = task.ext.args ?: ""
    """
    samtools mpileup ${bam} $args  > all.tsv
    samtools mpileup ${bam1} $args > term1.tsv
    samtools mpileup ${bam3} $args > term3.tsv
    """
}