process BAM_RMDUP {
    container (workflow.containerEngine ? "merszym/biohazard_bamrmdup:v0.2.2" : null)
    tag "${meta.RG}:${meta.Family}:${meta.Species}"
    label "process_low"
    label 'local'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("deduped_${bam}"), emit: bam
    tuple val(meta), path("rmdup.txt")     , emit: txt
    path "versions.yml"                    , emit: versions

    script:
    out_bam = "${meta.Taxon}/${meta.Reference}/3-deduped/${meta.RG}.${meta.Family}.${meta.Species}_deduped.bam"
    """
    bam-rmdup -r -o \"deduped_${bam}\" \"${bam}" > rmdup.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bam-rmdup --version 2>&1> /dev/null | cut -d ' ' -f3
    END_VERSIONS
    """
}