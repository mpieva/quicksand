#!/usr/bin/env nextflow

include { SPLITBAM }        from './modules/local/splitbam'
include { ESTIMATE_CC }     from './modules/local/ccestimate'

workflow{
    versions = Channel.empty()
    bam = file(params.bam)
    by  = file(params.by)

    SPLITBAM([['id':'test'],bam,by])
    versions = versions.mix(SPLITBAM.out.versions)

    ESTIMATE_CC(SPLITBAM.out.stats)

    ESTIMATE_CC.out.txt
        .filter{it[1].text != ''}
        .collectFile(storeDir: '.', newLine:true) {meta, file ->
            [ "cc_estimates.txt", ["# Cross contamination calculated from File: ${meta.RG}.txt"  ,file.text].join('\n') ]
        }

}
