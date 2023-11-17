include { SAMTOOLS_FASTA       } from '../modules/local/samtools_fasta'
include { RUN_KRAKENUNIQ       } from '../modules/local/krakenuniq_run'
include { PARSE_KRAKEN_REPORT  } from '../modules/local/parse_kraken_report'

workflow krakenrun {
    take: bam
    take: database

    main:
        // Convert bam to fasta
        SAMTOOLS_FASTA( bam )
        versions = SAMTOOLS_FASTA.out.versions.first()

        // Run krakenuniq
        RUN_KRAKENUNIQ( SAMTOOLS_FASTA.out.fasta.combine(database) )
        versions = versions.mix( RUN_KRAKENUNIQ.out.versions.first() )

        // Parse the krakenuniq-report
        PARSE_KRAKEN_REPORT( RUN_KRAKENUNIQ.out.report )
        versions = versions.mix( PARSE_KRAKEN_REPORT.out.versions.first() )

        // prepare the output channel for best_species
        PARSE_KRAKEN_REPORT.out.parsed_report
            .map{ meta, report ->
                [meta, report.splitCsv(sep:'\t', header:true)]
            }
            .branch{
                empty: it[1].size() == 0
                assigned: true
            }
            .set{parsed_report}

        // Handle Readgroups with no assignment
        parsed_report.empty.set{ empty }

        // Handle Readgroups with assignment(s)
        // Include a key [RG, taxon] to merge references later back to the extracted bams
        parsed_report.assigned
            .transpose()
            .map{meta, report ->
                [[meta.id, report[['f':'Family','o':'Order'][params.taxlvl]]], meta, report] // extract the 'taxon' from the parsed report
            }
            .set{assignments}

        // parse and filter the translate-file
        translate = RUN_KRAKENUNIQ.out.translate
        translate.map{ meta, translate ->
            [meta, translate, translate.readLines()]
        }
        .transpose()
        .filter{ it[2] =~ "${params.taxlvl}__" }
        .map{ meta, translate, asgn ->
            [meta + ['Taxon':(asgn =~ "${params.taxlvl}__([^|]*)")[0][1]],
            translate]
        }
        .unique()
        .map{ meta, translate ->
            [[meta.id, meta.Taxon], meta, translate ]
        }
        .combine(assignments, by:0) // this is to filter for the taxa in the parsed_report
        .map{ key, meta, translate, meta2, report ->
            [ meta, translate ]
        }
        .set{ translate }

    emit:
        versions = versions
        translate = translate
        assignments = assignments
        empty = empty
}