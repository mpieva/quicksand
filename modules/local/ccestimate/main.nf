process ESTIMATE_CC{
    container (workflow.containerEngine ? "pypy:3" : null)
    label 'local'

    input:
    tuple val(meta), path(stats)

    output:
    path "versions.yml"            , emit: versions
    path "cc_estimates.txt"

    script:
    """
    cross_cont.py ${stats} > cc_estimates.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypy: \$(pypy3 --version | tail -1 | cut -d ' ' -f2)
    END_VERSIONS
    """
}