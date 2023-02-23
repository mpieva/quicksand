#!/usr/bin/env nextflow

// include workflows for different executions of the pipeline
include { splitbam }   from './workflows/splitbam'
include { splitdir }   from './workflows/splitdir'
include { bamfilter }  from './workflows/bamfilter'

// include modules that are used by the main workflow

include { SAMTOOLS_FASTA } from './modules/local/samtools_fasta'
include { RUN_KRAKENUNIQ } from './modules/local/krakenuniq_run'


// input
versions = Channel.empty()
bam   = params.bam   ? file( params.bam,   checkIfExists:true) : ""
by    = params.by    ? file( params.by,    checkIfExists:true) : ""
split = params.split ? Channel.fromPath("${params.split}/*", checkIfExists:true) : ""

//
//
// The main workflow
//
//

workflow {

    //
    // Input Processing ~ Input Parameters
    //

    if (bam) {
        splitbam( bam,by )

        bam = splitbam.out.bams
        cc_stats = splitbam.out.cc
        versions = versions.mix( splitbam.out.versions )
    }
    else {
        splitdir( split )

        bam = splitdir.out.bams
        cc_stats = splitdir.out.cc
        versions = versions.mix( splitdir.out.versions )
    }

    //
    // Save the crosscontamination file
    //

    cc_stats
        .map{ it[1].text }
        .collectFile( storeDir:'.', newLine:true, name:"cc_estimates.txt" )

    //
    // Filter the bam files
    //

    bam.map { [it[0] + [ "id":it[1].baseName], it[1]] }.set{ bam }
    bamfilter( bam )

    bam = bamfilter.out.bam
    versions = versions.mix( bamfilter.out.versions )

    //
    // Convert bam to fasta for kraken
    //

    SAMTOOLS_FASTA( bam )
    versions = versions.mix( SAMTOOLS_FASTA.out.versions.first() )

    //
    // Run kraken
    //

    database = Channel.fromPath("${params.db}", type:'dir', checkIfExists:true)
    for_kraken = SAMTOOLS_FASTA.out.fasta.combine(database)

    RUN_KRAKENUNIQ( for_kraken )
    versions = versions.mix( RUN_KRAKENUNIQ.out.versions.first() )

}
