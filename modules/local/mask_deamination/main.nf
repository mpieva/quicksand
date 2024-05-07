process MASK_DEAMINATION{
    container (workflow.containerEngine ? "merszym/bam_deam:nextflow" : null)
    tag "${meta.id}:${meta.Taxon}:${meta.Species}"
    label 'local'

    input:
    tuple val(meta), path(bam1), path(bam3), path(bam)

    output:
    tuple val(meta), path("masked_${bam1}"), path("masked_${bam3}"), path("masked_${bam}"), emit: bam

    script:
    """
    mask_qual_scores.py ${bam} ${params.doublestranded ? 'doublestranded' : ''}
    mask_qual_scores.py ${bam1} ${params.doublestranded ? 'doublestranded' : ''}
    mask_qual_scores.py ${bam3} ${params.doublestranded ? 'doublestranded' : ''}
    """
}
