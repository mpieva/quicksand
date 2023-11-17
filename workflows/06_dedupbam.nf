include { BAM_RMDUP         } from '../modules/local/bam_rmdup'
include { SAMTOOLS_COVERAGE } from '../modules/local/samtools_coverage'

workflow dedupbam {
    take: mapped_bam
    main:
        // Take the mapped bam files
        // remove deduplicates

        BAM_RMDUP( mapped_bam )

        // Add the Number of Unique Sequences to the meta
        // first parse the bam-rmdup output file
        stats = BAM_RMDUP.out.txt.splitCsv(sep:'\t', header:true, limit:1)

        //now add the stats to the meta
        BAM_RMDUP.out.bam
        .combine(stats, by: 0)
        .map{ meta, bam, stats ->
            [
                meta+["ReadsDeduped":stats["out"].replace(",","") as int],
                bam
            ]
        }
        .set{ bam }

        versions = BAM_RMDUP.out.versions.first()

        // Calculate the coverage

        SAMTOOLS_COVERAGE( bam )
        versions = versions.mix(SAMTOOLS_COVERAGE.out.versions.first())

        // Add the CoveredBP value to the meta of the main bam channel

        SAMTOOLS_COVERAGE.out.bam
        .map{ meta, bam, cov ->
            [
                meta+["CoveredBP":cov.trim() as int],
                bam
            ]
        }
        .set{ bam }

    emit:
        bam = bam
        versions = versions
}