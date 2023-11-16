include { BEDTOOLS_INTERSECT } from '../modules/local/bedtools_intersect'
include { SAMTOOLS_COUNT     } from '../modules/local/samtools_count'
include { SAMTOOLS_COVERAGE  } from '../modules/local/samtools_coverage'

workflow bedfilterbam {
    take: best
    take: bedfiles
    main:
        // Take the "best subset"
        // find the species bed file in the masked folder

        bedfiles.map{ bedfile ->
            def matcher = bedfile.name =~ /(.+)\.masked\.bed/
            def speciesName = matcher ? matcher[0][1] : null
            [speciesName, bedfile]
        }.set{ bedfiles }

        // now combine with the bedfile
        // so that bedtools intersect works
        // first, create the species as key
        best.map{ meta, bam ->
            [ meta.Species, meta, bam ]
        }
        .combine( bedfiles, by:0 ) // [sp, meta, bam, bed] --> now remove the sp again
        .map{ sp, meta, bam, bed ->
            [ meta, bam, bed ]
        }
        .set{ best }

        // And run bedtools instersect
        BEDTOOLS_INTERSECT( best )

        best = BEDTOOLS_INTERSECT.out.bam
        versions = BEDTOOLS_INTERSECT.out.versions.first()

        // get the counts
        SAMTOOLS_COUNT( best )
        // add the counts to the meta
        best = SAMTOOLS_COUNT.out.bam.map { meta, bam, count ->
            [ meta+[ 'ReadsBedfiltered': count as int ], bam ]
        }

        // get the coverage
        SAMTOOLS_COVERAGE( best )
        // And add the coverage to the meta
        SAMTOOLS_COVERAGE.out.bam
        .map{ meta, bam, cov ->
            [
                meta+[ "PostBedCoveredBP":cov.trim() as int ],
                bam
            ]
        }
        .set{ best }

    emit:
        bam = best
        versions = versions
}