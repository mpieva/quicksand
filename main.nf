#!/usr/bin/env nextflow

// include workflows for different executions of the pipeline
include { setup      } from './workflows/00_setup'
include { splitbam   } from './workflows/01_splitbam'
include { splitdir   } from './workflows/01_splitdir'
include { bamfilter  } from './workflows/02_bamfilter'
include { bamextract } from './workflows/04_bamextract'
include { krakenrun  } from './workflows/03_krakenrun'
include { refprep    } from './workflows/03_refprep'
include { mapbam     } from './workflows/05_mapbam'
include { dedupbam   } from './workflows/06_dedupbam'

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
    // 0. Setup the folders etc.
    //

    setup([])

    //
    // 1. Input Processing ~ Input Parameters
    //

    if (bam) {
        splitbam( bam,by )

        bam = splitbam.out.bams
        versions = versions.mix( splitbam.out.versions )
    }
    else {
        splitdir( split )

        bam = splitdir.out.bams
        versions = versions.mix( splitdir.out.versions )
    }

    //
    // 2. Filter the bam files
    //

    bam.map { [it[0] + [ "id":it[1].baseName, 'Reference':'best'], it[1]] }.set{ bam }
    bamfilter( bam )

    bam = bamfilter.out.bam
    versions = versions.mix( bamfilter.out.versions )

    //
    // 3. Run kraken
    //

    krakenrun( bam, database )
    // #TODO: handle empty...
    version = versions.mix( krakenrun.out.versions )

    //
    // 4.1 Extract bams based on kraken-results
    //

    bamextract( bamfilter.out.bam, krakenrun.out.translate )
    versions = versions.mix( bamextract.out.versions )

    //
    // 4.2 Prepare the reference genomes
    //

    assignments = krakenrun.out.assignments

    refprep( database, assignments, [] )
    versions = versions.mix( refprep.out.versions )

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
    .set{bwa_in}

    //
    // 5. Map with BWA
    //

    mapbam( bwa_in )
    versions = versions.mix( mapbam.out.versions )

    //
    // 6. Dedup the mapped bam
    //

    dedupbam(mapbam.out.bam)
    versions = versions.mix( dedupbam.out.versions )
}
