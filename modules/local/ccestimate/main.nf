process ESTIMATE_CC{
    container (workflow.containerEngine ? "merszym/bam_deam:nextflow" : null)
    label 'local'

    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path("CC.txt"), emit: txt

    script:
    """
    cross_cont.py ${stats} > CC.txt
    """
}