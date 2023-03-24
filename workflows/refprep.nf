include { PARSE_TAXONOMY } from '../modules/local/parse_taxonomy'

workflow get_reference {
    take: database
    take: assignments
    take: fixed_references
    take: genomes

    main:
    PARSE_TAXONOMY( database )
    json = PARSE_TAXONOMY.out.json

}