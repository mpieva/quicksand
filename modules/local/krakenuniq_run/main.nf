process RUN_KRAKENUNIQ {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakenuniq:0.7.3--pl5321h19e8d03_0' :
        'quay.io/biocontainers/krakenuniq:0.7.3--pl5321h19e8d03_0' }"

    tag "$meta.id"

    input:
    tuple val(meta), path(fasta), path(database)

    output:
    tuple val(meta), path("${meta.id}.kraken.translate"), emit: translate
    tuple val(meta), path("${meta.id}.kraken.report")   , emit: report
    path "versions.yml"                                 , emit: versions

    script:
    """
    krakenuniq --threads ${task.cpus} \
    --db ${database} \
    --fasta-input ${fasta} \
    --report-file ${meta.id}.kraken.report > output.kraken

    krakenuniq-translate \
    --db ${database} \
    --mpa-format output.kraken > ${meta.id}.kraken.translate

    echo "\$(grep -E '^(#|%|[0-9]).*' ${meta.id}.kraken.report)" > ${meta.id}.kraken.report

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakenuniq: \$(krakenuniq --version | head -1 | cut -f3 -d ' ')
    END_VERSIONS
    """
}