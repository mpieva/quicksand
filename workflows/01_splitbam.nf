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
        cc_versions = ESTIMATE_CC.out.versions

    emit:
        versions = cc_versions.mix(SPLITBAM.out.versions)
        bams = SPLITBAM.out.bams.transpose()
}
