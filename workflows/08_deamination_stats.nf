include { BAM_DEAM_STATS as BAM_DEAM_BEST  } from '../modules/local/bam_deam_stats'
include { BAM_DEAM_STATS as BAM_DEAM_FIXED } from '../modules/local/bam_deam_stats'
include { MASK_DEAMINATION                 } from '../modules/local/mask_deamination'
include { SAMTOOLS_MPILEUP                 } from '../modules/local/samtools_mpileup'
include { PLOT_DEAM                        } from '../modules/local/pandas_plot_deam'

workflow deamination_stats {
    take: best
    take: fixed
    main:
        // Take the "best subset"
        // Get the deamination stats
        // and be done with it

        BAM_DEAM_BEST( best )

        BAM_DEAM_BEST.out.tsv
        .map{ meta, tsv ->
            def deam = tsv.splitCsv(sep:'\t', header:true).first() // first because the splitCsv results have only one row -> [[key:value]]
            [
                meta+deam,
            ]
        }
        .set{ best }

        ch_positions = BAM_DEAM_BEST.out.positions

        // Now for the fixed

        BAM_DEAM_FIXED( fixed )
        // get the deaminated bam-files
        fixed = BAM_DEAM_FIXED.out.bam
        // and get the stats
        fixed.combine( BAM_DEAM_FIXED.out.tsv, by:0 )
        .map{ meta, bam1, bam3, bam, tsv ->
            def deam = tsv.splitCsv(sep:'\t', header:true).first() // first because the splitCsv results in [[key:value]]
            [
                meta+deam,
                bam1,
                bam3,
                bam
            ]
        }
        .set{ fixed }

        ch_positions = ch_positions.mix(BAM_DEAM_FIXED.out.positions)
        
        // Plot deamination rates
        PLOT_DEAM(ch_positions)

        // And continue with the fixed files...

        MASK_DEAMINATION( fixed )

        //
        // Create Mpileups
        //

        SAMTOOLS_MPILEUP( MASK_DEAMINATION.out.bam )


    emit:
        best = best
        fixed = SAMTOOLS_MPILEUP.out.tsv
}