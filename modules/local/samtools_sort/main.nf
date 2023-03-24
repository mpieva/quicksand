process SAMTOOLS_SORT{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    tag "$meta.id:$meta.taxon"
    label "process_low"
    label "local"
    publishDir 'out', mode: 'copy', saveAs: {out_bam}, pattern: 'sorted_*.bam'

    input:
    tuple val(meta), path(extracted_bam)

    output:
    tuple val(meta), path("sorted_${extracted_bam}"), emit: bam
    path "versions.yml"                             , emit: versions

    script:
    out_bam = "${meta.taxon}/1-extracted/${meta.id}_extractedReads-${meta.taxon}.bam"
    """
    samtools sort -n -l $params.compression_level -o sorted_${extracted_bam}  ${extracted_bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """
}