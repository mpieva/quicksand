// running in batch and preloading the database doubled the time for me in the tests??
// I keep the process here for now, but I go back to the old 1:1 kraken-runs
// doesnt make sense, but maybe its the nextflow overhead?

process RUN_KRAKENUNIQ_BATCH {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakenuniq:1.0.2--pl5321h19e8d03_0':
        'quay.io/biocontainers/krakenuniq:1.0.2--pl5321h19e8d03_0' }"

    input:
    path(fasta)
    path database

    output:
    path "*.kraken.translate", emit: translate
    path "*.kraken.report"   , emit: report
    path "versions.yml"      , emit: versions

    script:
    """
    krakenuniq --db ${database} \
    --threads ${task.cpus} \
    --preload

    printf "%s\\n" ${fasta} | while read FASTA; do \\
        PREFIX="\${FASTA%%.*}"

        krakenuniq --threads ${task.cpus} \
        --db ${database} \
        --report-file \${PREFIX}.kraken.report \
        --output output.kraken \
        \${FASTA}

        krakenuniq-translate \
        --db ${database} \
        --mpa-format output.kraken > \${PREFIX}.kraken.translate
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakenuniq: \$(krakenuniq --version | head -1 | cut -f3 -d ' ')
    END_VERSIONS
    """
}

process RUN_KRAKENUNIQ {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakenuniq:1.0.2--pl5321h19e8d03_0':
        'quay.io/biocontainers/krakenuniq:1.0.2--pl5321h19e8d03_0' }"

    input:
    tuple val(meta), path(fasta), path(database)

    output:
    path "${meta.id}.kraken.translate", emit: translate
    path "${meta.id}.kraken.report"   , emit: report
    path "versions.yml"               , emit: versions

    script:
    """
    krakenuniq --threads ${task.cpus} \
    --db ${database} \
    --report-file ${meta.id}.kraken.report \
    --output output.kraken \
    ${fasta}

    krakenuniq-translate \
    --db ${database} \
    --mpa-format output.kraken > ${meta.id}.kraken.translate

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakenuniq: \$(krakenuniq --version | head -1 | cut -f3 -d ' ')
    END_VERSIONS
    """
}