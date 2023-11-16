include { PARSE_TAXONOMY } from '../modules/local/parse_taxonomy'

/*
Prepare the reference-files for mapping with BWA
1. get for each taxon id the reference-genomes
2. replace the references if the references is fixed for a certain taxon
*/

import groovy.json.JsonSlurper

workflow refprep {
    take: database
    take: assignments
    take: fixed_references

    main:
    //
    // 1. Parse the taxonomy
    //

    def jsonSlurper = new JsonSlurper()

    PARSE_TAXONOMY( database )
    PARSE_TAXONOMY.out.versions.set{ versions }
    PARSE_TAXONOMY.out.json.map{ json ->
        [jsonSlurper.parseText(json.text)]
    }
    .set{ json }

    //
    // 2. Get for each 'extracted bam' the list of references
    //

    assignments.combine(json)
    .map{key, meta, report, json ->
        [key, report, json[report.BestTaxID]]
    }
    .set{ assignments }

    //
    // #TODO: 3. Handle fixed references
    //

    emit:
        references = assignments
        versions = versions

}