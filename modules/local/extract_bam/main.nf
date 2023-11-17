process EXTRACT_TAXON {
    container (workflow.containerEngine ? "merszym/bamfilter:nextflow" : null)
    tag "$meta.id:$meta.Taxon"
    label "process_low"
    label 'local'

    input:
    tuple val(meta), path(txt), path(bam)

    output:
    tuple val(meta), path("${meta.id}_${meta.Taxon}.out.bam"), emit: bam
    path "versions.yml"                                      , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    bamfilter -i ${txt} $args -o ${meta.id}_${meta.Taxon}.out.bam ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamfilter: \$(bamfilter --version | cut -f2 -d ' ')
    END_VERSIONS
    """
}