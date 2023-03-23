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

    emit:
        versions = versions
        translate = RUN_KRAKENUNIQ.out.translate
        report = PARSE_KRAKEN_REPORT.out.parsed_report
}