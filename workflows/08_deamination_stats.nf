include { BAM_DEAM_STATS as BAM_DEAM_BEST  } from '../modules/local/bam_deam_stats'
include { BAM_DEAM_STATS as BAM_DEAM_FIXED } from '../modules/local/bam_deam_stats'
include { SAMTOOLS_COUNT                   } from '../modules/local/samtools_count'
include { SAMTOOLS_COVERAGE                } from '../modules/local/samtools_coverage'

workflow deamination_stats {
    take: best
    take: fixed
    main:
        // Take the "best subset"
        // find the species bed file in the masked folder

        BAM_DEAM_BEST( best )

        best.combine( BAM_DEAM_BEST.out.tsv, by:0 )
        .map{ meta, bam, tsv ->
            def deam = tsv.splitCsv(sep:'\t', header:true).first() // first because the splitCsv results in [[key:value]]
            [
                meta+deam,
                bam
            ]
        }
        .set{ best }

        // Now for the fixed

        BAM_DEAM_FIXED( fixed )
        // get the deaminated bam-files
        fixed = BAM_DEAM_FIXED.out.bam
        // and get the stats
        fixed.combine( BAM_DEAM_FIXED.out.tsv, by:0 )
        .map{ meta, bam1, bam3, tsv ->
            def deam = tsv.splitCsv(sep:'\t', header:true).first() // first because the splitCsv results in [[key:value]]
            [
                meta+deam,
                bam1,
                bam3
            ]
        }
        .set{ fixed }

    emit:
        best = best
}