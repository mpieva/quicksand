process ESTIMATE_CC{
    container (workflow.containerEngine ? "merszym/bam_deam:nextflow" : null)
    label 'local'

    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path("${meta.RG}.CC.txt"), emit: txt

    script:
    """
    cross_cont.py splittingstats.txt > \"${meta.RG}.CC.txt\"
    """
}