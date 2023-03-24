process EXTRACT_TAXON {
    container (workflow.containerEngine ? "merszym/bamfilter:nextflow" : null)
    tag "$meta.id:$meta.taxon"
    label "process_low"
    label 'local'

    input:
    tuple val(meta), path(txt), path(bam)

    output:
    tuple val(meta), path("${meta.id}_${meta.taxon}.out.bam"), emit: bam
    path "versions.yml"                                      , emit: versions

    when:
    meta.ReadsExtracted > params.krakenuniq_min_reads

    script:
    """
    bamfilter -i ${txt} -l ${params.compression_level} -o ${meta.id}_${meta.taxon}.out.bam ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamfilter: \$(bamfilter --version | cut -f2 -d ' ')
    END_VERSIONS
    """
}