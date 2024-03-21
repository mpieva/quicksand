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

        // Add the samtools coverage values to the meta of the main bam channel

        SAMTOOLS_COVERAGE.out.bam
        .map{ meta, bam, cov ->
            def (covered_bases, b, c) = cov.split()
            def coverage = c.trim() as float
            def breadth = (b.trim() as float) / 100
            def expected_breadth = 1 - (Math.exp(-0.883 * coverage))
            [
                meta+[
                        "CoveredBP": covered_bases.trim() as int,
                        "Coverage": coverage.trunc(2),
                        "Breadth": breadth.trunc(3),
                        'ExpectedBreadth': expected_breadth.trunc(3),
                        'ProportionExpectedBreadth': (breadth / expected_breadth).trunc(3)
                    ],
                bam
            ]
        }
        .set{ bam }

    emit:
        bam = bam
        versions = versions
}