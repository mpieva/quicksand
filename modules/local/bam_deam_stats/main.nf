process BAM_DEAM_STATS{
    container (workflow.containerEngine ? "merszym/bam_deam:nextflow" : null)
    tag "${meta.id}:${meta.Taxon}:${meta.Species}"
    label 'local'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('ancient_stats.tsv')                                                , emit: tsv
    tuple val(meta), path("output.deaminated1.bam"), path("output.deaminated3.bam"), path(bam), emit: bam

    script:
    def args = task.ext.args ?: ''
    """
    touch output.deaminated1.bam
    touch output.deaminated3.bam

    bam_deam_stats.py ${bam} $args > ancient_stats.tsv
    """
}