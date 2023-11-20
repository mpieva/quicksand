process BAM_RMDUP {
    container (workflow.containerEngine ? "merszym/biohazard_bamrmdup:v0.2.2" : null)
    tag "${meta.id}:${meta.Family}:${meta.Species}"
    label "process_low"
    label 'local'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("deduped_${bam}"), emit: bam
    tuple val(meta), path("rmdup.txt")     , emit: txt
    path "versions.yml"                    , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    bam-rmdup $args -o deduped_${bam} ${bam} > rmdup.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bam-rmdup: \$(bam-rmdup --version 2>&1> /dev/null | cut -d ' ' -f3)
    END_VERSIONS
    """
}