process ESTIMATE_CC{
    container (workflow.containerEngine ? "pypy:3" : null)
    label 'local'

    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path("CC.txt"), emit: txt
    path "versions.yml"            , emit: versions

    script:
    """
    cross_cont.py ${stats} > CC.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypy: \$(pypy3 --version | tail -1 | cut -d ' ' -f2)
    END_VERSIONS
    """
}