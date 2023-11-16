process MAP_BWA {
    container (workflow.containerEngine ? "merszym/network-aware-bwa:v0.5.10" : null)
    tag "${meta.id}:${meta.Taxon}:${meta.Species}"
    label "process_low"
    label 'local'

    input:
    tuple val(meta), path(bam), path(genomesdir)

    output:
    tuple val(meta), path("mapped_${bam}"), emit: bam
    path 'versions.yml'                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    mapbwa_out = "${meta.Taxon}/${meta.Reference}/2-aligned/${meta.id}.${meta.Family}.${meta.Species}.bam"
    println mapbwa_out
    """
    bwa bam2bam -g \"${genomesdir}/${meta.Family}/${meta.Species}.fasta\" $args --only-aligned ${bam} > mapped_${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 > /dev/null | grep Version | cut -d ' ' -f2)
    END_VERSIONS
    """
}