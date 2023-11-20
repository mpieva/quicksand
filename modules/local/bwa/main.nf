process MAP_BWA {
    container (workflow.containerEngine ? "merszym/network-aware-bwa:v0.5.10" : null)
    tag "${meta.id}:${meta.Taxon}:${meta.Species}"
    label "process_low"
    label 'local'

    input:
    tuple val(meta), path(bam), path(genome)

    output:
    tuple val(meta), path("mapped_${bam}"), emit: bam
    path 'versions.yml'                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    bwa index ${genome}
    bwa bam2bam -g ${genome} $args --only-aligned ${bam} > mapped_${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 > /dev/null | grep Version | cut -d ' ' -f2)
    END_VERSIONS
    """
}