process PLOT_DEAM{
    container (workflow.containerEngine ? "merszym/pandas_plt:nextflow" : null)
    tag "${meta.id}"
    label 'local'

    input:
    tuple val(meta), path(tsv), path(confidence)

    output:
    path("*.jpg") , emit: jpg

    script:
    """
    plot_deam_patterns.py ${tsv} ${meta.Family}.${meta.Species}
    """
}