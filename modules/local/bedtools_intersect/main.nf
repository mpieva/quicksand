process BEDTOOLS_INTERSECT {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3' :
        'quay.io/biocontainers/bedtools:2.30.0--h468198e_3' }"
    label "process_low"
    label 'local'
    tag "${meta.id}:${meta.Family}:${meta.Species}"

    input:
    tuple val(meta), path(bam), path(bedfile)

    output:
    tuple val(meta), path("masked_${bam}"), emit: bam
    path "versions.yml"                   , emit: versions

    script:
    """
    bedtools intersect -a ${bam} -b ${bedfile} -v > masked_${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | cut -d' ' -f2)
    END_VERSIONS
    """
}