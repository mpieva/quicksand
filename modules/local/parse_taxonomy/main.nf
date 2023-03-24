process PARSE_TAXONOMY{
    container (workflow.containerEngine ? "pypy:3" : null)
    label '$database'
    label 'local'

    input:
    path(database)

    output:
    path("taxonomy.json"), emit: json

    script:
    """
    parse_taxonomy.py ${database}/taxonomy/nodes.dmp ${database}/taxonomy/names.dmp ${database}/seqid2taxid.map > taxonomy.json
    """
}