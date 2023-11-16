process BAM_LENGTHFILTER {
    container (workflow.containerEngine ? "merszym/bam-lengthfilter:nextflow" : null)
    tag "$meta.id"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("length${bam}"), emit: bam
    path "versions.yml"                  , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    bam-lengthfilter \
      $args \
      -o length${bam} \
      ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bam-lengthfilter: \$(bam-lengthfilter --version | cut -d' ' -f2)
    END_VERSIONS
    """
}