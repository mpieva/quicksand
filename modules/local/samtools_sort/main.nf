process SAMTOOLS_SORT{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    tag "$meta.id:$meta.Taxon"
    label "local"

    input:
    tuple val(meta), path(extracted_bam)

    output:
    tuple val(meta), path("sorted_${extracted_bam}"), emit: bam
    path "versions.yml"                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    samtools sort $args -o sorted_${extracted_bam}  ${extracted_bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """
}