#!/usr/bin/env nextflow

// include workflows for different executions of the pipeline
include { splitbam }  from './workflows/splitbam'
include { splitdir }  from './workflows/splitdir'

// include modules that are used by the main workflow
include { SAMTOOLS_FILTER } from './modules/local/samtools_filter'

// some required functions
def has_ending(file, extension){
    return extension.any{ file.toString().toLowerCase().endsWith(it) }
}

// input
versions = Channel.empty()
bam   = params.bam   ? file( params.bam,   checkIfExists:true) : ""
by    = params.by    ? file( params.by,    checkIfExists:true) : ""
split = params.split ? Channel.fromPath("${params.split}/*", checkIfExists:true) : ""


// The main workflow

workflow {
    if (bam) {
        splitbam( bam,by )

        bams = splitbam.out.bams
        cc_stats = splitbam.out.cc
        versions = versions.mix( splitbam.out.versions )
    }
    else {
        splitdir( split )

        bams = splitdir.out.bams
        cc_stats = splitdir.out.cc
    }
    // save the crosscontamination file
    cc_stats
        .filter { it[1].text != '' }
        .collectFile( storeDir: '.', newLine:true ) { meta, file ->
            [ "cc_estimates.txt", ["# Cross contamination calculated from File: ${meta.RG}.txt"  ,file.text ].join( '\n' ) ]
        }

    bams.map { [it[0] + [ "id":it[1].baseName], it[1]] }.set{ bams }

    // filter bam
    SAMTOOLS_FILTER(bams)
}
