process CHECK_DATABASE {
    container (workflow.containerEngine ? "pypy:3" : null)
    label 'local'
    label "process_low"

    input:
    path(database)
    path(genomes)
    path(bedfiles)

    output:
    path("database_report.tsv"), emit: report
    path("summary.txt")        , emit: summary

    script:
    """
    check_database.sh -d ${database} -g ${genomes} -b ${bedfiles} > database_report.tsv 2> summary.txt
    """
}
