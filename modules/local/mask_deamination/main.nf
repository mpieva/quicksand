process MASK_DEAMINATION{
    container (workflow.containerEngine ? "merszym/bam_deam:nextflow" : null)
    tag "${meta.id}:${meta.Taxon}:${meta.Species}"
    label 'process_medium'
    label 'local'

    input:
    tuple val(meta), path(bam1), path(bam3), path(bam)

    output:
    tuple val(meta), path("masked_${bam}"), path("masked_${bam1}"), path("masked_${bam3}"), emit: bam

    script:
    """
    mask_qual_scores.py ${bam}
    mask_qual_scores.py ${bam1}
    mask_qual_scores.py ${bam3}
    """
}
