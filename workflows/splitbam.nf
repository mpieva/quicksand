#!/usr/bin/env nextflow

include { SPLITBAM }        from '../modules/local/splitbam'
include { ESTIMATE_CC }     from '../modules/local/ccestimate'

workflow splitbam {
    take:
        bam
        by
    main:
        // Split Bam by RGs
        SPLITBAM( [[:], bam, by] )

        // Estimate Crosscontamination
        ESTIMATE_CC(SPLITBAM.out.stats)

    emit:
        versions = SPLITBAM.out.versions
        bams = SPLITBAM.out.bams.transpose()
        cc = ESTIMATE_CC.out.txt
}
