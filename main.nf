#!/usr/bin/env nextflow

include { SPLITBAM }        from './modules/local/splitbam'
include { ESTIMATE_CC }     from './modules/local/ccestimate'
include { SAMTOOLS_FQ2BAM } from './modules/local/samtools_fq2bam'


// some required functions
def has_ending(file, extension){
    return extension.any{ file.toString().toLowerCase().endsWith(it) }
}

// input
versions = Channel.empty()
bam   = params.bam   ? file( params.bam,   checkIfExists:true) : ""
by    = params.by    ? file( params.by,    checkIfExists:true) : ""
split = params.split ? Channel.fromPath("${params.split}/*", checkIfExists:true) : ""


// Define all the alternative workflows

workflow splitdir {
    take:
        split
    main:
        split
        .map{ [[:], it] }
            .branch {
                bam: it[1].getExtension() == 'bam'
                fastq: has_ending( it[1], ["fastq","fastq.gz","fq","fq.gz"])
                stats:  it[1].name =~ /split.*stats/
                fail: true
            }
            .set{ split }

        split.fail.view()
        // Estimate cross-contamination if file exists
        ESTIMATE_CC( split.stats.first() )

        // convert fastq to bam
        SAMTOOLS_FQ2BAM( split.fastq )

    emit:
        bams = split.bam.mix( SAMTOOLS_FQ2BAM.out.bam )
        cc = ESTIMATE_CC.out.txt
}


workflow splitbam {
    take:
        bam
        by
    main:
        // Split Bam by RGs
        SPLITBAM( [[:], bam, by] )
        versions = versions.mix( SPLITBAM.out.versions )

        // Estimate Crosscontamination
        ESTIMATE_CC(SPLITBAM.out.stats)

    emit:
        bams = SPLITBAM.out.bams.transpose()
        cc = ESTIMATE_CC.out.txt
}

// The main workflow

workflow {
    if (bam) {
        splitbam( bam,by )

        bams = splitbam.out.bams
        cc_stats = splitbam.out.cc
    }
    else {
        splitdir( split )

        bams = splitdir.out.bams
        cc_stats = splitdir.out.cc
    }
    cc_stats
        .filter { it[1].text != '' }
        .collectFile( storeDir: '.', newLine:true ) { meta, file ->
            [ "cc_estimates.txt", ["# Cross contamination calculated from File: ${meta.RG}.txt"  ,file.text ].join( '\n' ) ]
        }
}
