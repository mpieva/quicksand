#!/usr/bin/env nextflow

// include workflows for different executions of the pipeline
include { splitbam   } from './workflows/splitbam'
include { splitdir   } from './workflows/splitdir'
include { bamfilter  } from './workflows/bamfilter'
include { bamextract } from './workflows/bamextract'
include { krakenrun  } from './workflows/krakenrun'
include { refprep    } from './workflows/refprep'
include { mapbam     } from './workflows/mapbam'

// input
versions = Channel.empty()
bam        = params.bam     ? file( params.bam, checkIfExists:true) : ""
by         = params.by      ? file( params.by,  checkIfExists:true) : ""
split      = params.split   ? Channel.fromPath("${params.split}/*", checkIfExists:true) : ""
genomesdir = params.genomes ? Channel.fromPath("${params.genomes}", type:'dir',   checkIfExists:true)  : Channel.empty()


database = Channel.fromPath("${params.db}", type:'dir', checkIfExists:true)
taxid = new File("${params.genomes}/taxid_map.tsv").exists() ? Channel.fromPath("${params.genomes}/taxid_map.tsv", type:'file') : Channel.fromPath("${baseDir}/assets/taxid_map_example.tsv", type:'file')

//
//
// The main workflow
//
//

workflow {

    //
    // 1. Input Processing ~ Input Parameters
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
    // 2. Save the crosscontamination file
    //

    cc_stats
        .map{ it[1].text }
        .collectFile( storeDir:'.', newLine:true, name:"cc_estimates.txt" )

    //
    // 3. Filter the bam files
    //

    bam.map { [it[0] + [ "id":it[1].baseName, 'Reference':'fixed'], it[1]] }.set{ bam }
    bamfilter( bam )

    bam = bamfilter.out.bam
    versions = versions.mix( bamfilter.out.versions )

    //
    // 4. Run kraken
    //

    krakenrun( bam, database )
    version = versions.mix( krakenrun.out.versions )

    //
    // 5.1 Extract bams based on kraken-results
    //

    bamextract( bamfilter.out.bam, krakenrun.out.translate )
    versions = versions.mix( bamextract.out.versions.first() )

    //
    // 5.2 Prepare the reference genomes
    //

    assignments = krakenrun.out.assignments

    refprep( database, assignments, [] )
    versions = versions.mix( refprep.out.versions.first() )

    // combine the extracted and assigned paths

    bamextract.out.bam.map{ meta, bam ->
        [[meta.id, meta.taxon], meta, bam]
    }
    .join( refprep.out.references )
    .map{ key, meta, bam, report, references ->
        [meta+report, bam, references]
    }
    .transpose()
    .map{ meta, bam, reference ->
        [meta+['Species':reference], bam]
    }
    .combine( genomesdir )
    .view()
    .set{bwa_in}

    //
    // 6. Map with BWA
    //

    mapbam( bwa_in )
    versions = versions.mix( mapbam.out.versions.first() )

}
