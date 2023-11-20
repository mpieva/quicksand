#!/usr/bin/env nextflow

// include workflows for different executions of the pipeline
include { setup             } from './workflows/00_setup'
include { splitbam          } from './workflows/01_splitbam'
include { splitdir          } from './workflows/01_splitdir'
include { bamfilter         } from './workflows/02_bamfilter'
include { bamextract        } from './workflows/04_bamextract'
include { krakenrun         } from './workflows/03_krakenrun'
include { refprep           } from './workflows/03_refprep'
include { mapbam            } from './workflows/05_mapbam'
include { dedupbam          } from './workflows/06_dedupbam'
include { bedfilterbam      } from './workflows/07_bedfilterbam'
include { deamination_stats } from './workflows/08_deamination_stats'
include { write_reports     } from './workflows/09_reports.nf'

//The colors
red = "\033[0;31m"
white = "\033[0m"
yellow = "\033[0;33m"

// Define some functions

def exit_with_error_msg(error, text){
    println "[quicksand]: ${red}${error}: ${text}${white}"
    exit 0
}
def get_warn_msg(text){
    return "[quicksand]: ${yellow}(WARN): ${text}${white}"
}
def get_info_msg(text){
    return "[quicksand]: ${text}"
}
def exit_missing_required(flag){
    exit_with_error_msg("ArgumentError", "missing required argument ${flag}")
}

// Define the current workflow

standard_run = params.rerun ? false : true

//
//
// Validation of input parameters
//
//

if(params.taxlvl !in ['f','o']){
    exit_with_error_msg("ArgumentError","taxlvl must be one of [o, f] not ${params.taxlvl}")
}
if( params.rerun && new File('final_report.tsv').exists()==false){
    log.info get_info_msg("Use --rerun only for additional runs together with --fixed flag")
    exit_with_error_msg("ArgumentError", "Conflicting arguments")
}
if(params.split && (params.bam || params.rg) && standard_run){
    log.info get_info_msg("Use: nextflow run mpieva/quicksand {--rg FILE --bam FILE | --split DIR}")
    exit_with_error_msg("ArgumentError", "Too many arguments")
}
if(!params.split && !(params.bam && params.rg) && standard_run){
    log.info get_info_msg("Use: nextflow run mpieva/quicksand {--rg FILE --bam FILE | --split DIR}")
    exit_with_error_msg("ArgumentError", "Too few arguments")
}
if(params.rerun && !(params.fixed)){
    log.info get_info_msg("Use --rerun together with --fixed")
    exit_with_error_msg("ArgumentError", "Too few arguments")
}
if(standard_run && !params.genomes){ exit_missing_required('--genomes') }
if(standard_run && !params.db){   exit_missing_required('--db')      }
if(standard_run && !params.bedfiles){ exit_missing_required('--bedfiles')}

//
//
// input Channels
//
//

bam        = params.bam      ? file( params.bam, checkIfExists:true) : ""
by         = params.rg       ? file( params.rg,  checkIfExists:true) : ""
split      = params.split    ? Channel.fromPath("${params.split}/*",     checkIfExists:true) : ""
genomesdir = params.genomes  ? Channel.fromPath("${params.genomes}/*/*.fasta", checkIfExists:true) : Channel.empty()
bedfiles   = params.bedfiles ? Channel.fromPath("${params.bedfiles}/*",  checkIfExists:true) : Channel.empty()

database = Channel.fromPath("${params.db}", type:'dir', checkIfExists:true)

// fixed references
ch_fixed = params.fixed ?
    Channel.fromPath("${params.fixed}", checkIfExists:true)
        .splitCsv(sep:'\t', header:['Family','Species','Genome'], skip:1)
    : Channel.empty()

ch_versions = Channel.empty()
ch_final = Channel.empty()

//
//
// The main workflow
//
//

workflow {

    standard_run = true

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
        ch_versions = ch_versions.mix( splitbam.out.versions )
    }
    else {
        splitdir( split )

        bam = splitdir.out.bams
        ch_versions = ch_versions.mix( splitdir.out.versions )
    }

    //
    // 2. Filter the bam files
    //

    //include a meta-file with all fields existing
    meta = Channel.fromPath("$baseDir/assets/pipeline/meta.tsv").splitCsv(sep:'\t', header:true)
    bam.combine(meta).map{ m1, bam, meta -> [meta, bam] }.set{ bam }

    bam.map {
        [
            it[0] + [
                "id":it[1].baseName,
                "RG":it[1].baseName,
                "ExtractLVL": params.taxlvl
            ],
            it[1]
        ]
    }.set{ bam }
    bamfilter( bam )

    bam = bamfilter.out.bam
    ch_versions = ch_versions.mix( bamfilter.out.versions )

    //
    // 3. Run kraken
    //

    krakenrun( bam, database )
    ch_versions = ch_versions.mix( krakenrun.out.versions )

    // Add the libraries with no assignments to the final channel
    ch_empty = krakenrun.out.empty.map{ it[0] }
    ch_final.mix(ch_empty).set{ ch_final }


    //
    // 4.1 Extract bams based on kraken-results
    //

    bamextract( bamfilter.out.bam, krakenrun.out.translate )
    ch_versions = ch_versions.mix( bamextract.out.versions )

    //
    // 4.2 Prepare the reference genomes
    //

    assignments = krakenrun.out.assignments

    refprep( database, assignments, [] )
    ch_versions = ch_versions.mix( refprep.out.versions )

    // combine the extracted and assigned paths

    bamextract.out.bam.map{ meta, bam ->
        [[meta.id, meta.Taxon], meta, bam]
    }
    .join( refprep.out.references )
    .map{ key, meta, bam, report, references ->
        [meta+report, bam, references]
    }
    .transpose()
    .map{ meta, bam, reference ->
        [meta+['Species':reference], bam]
    }
    .set{bwa_in}

    //
    // 5. Map with BWA
    //

    mapbam( bwa_in, genomesdir, ch_fixed )
    ch_versions = ch_versions.mix( mapbam.out.versions )

    //
    // 6. Dedup the mapped bam
    //

    dedupbam(mapbam.out.bam)
    ch_versions = ch_versions.mix( dedupbam.out.versions )

    deduped = dedupbam.out.bam

    // split between 'fixed' and 'best'

    deduped.branch{
        best:it[0].Reference == 'best'
        fixed:true
    }
    .set{deduped}

    // if default.best is empty it would throw an index error,
    // best: reduce the 1 "best" hit per family

    best = deduped.best
        .map{meta,bam -> [meta.id, meta.Family, meta.CoveredBP, meta, bam]}
        .toSortedList({ a,b -> a[0]+a[1] <=> b[0]+b[1] ?: a[2] <=> b[2]})
        .flatten()
        .collate(5)
        .groupTuple(by:[0,1])   //[[rg, fam, [covered_bp < .. < covered_bp][meta,meta,meta],[bam,bam,bam]]
        .map{n -> [n[3][-1], n[4][-1]]} // from the highest, the [meta, bam]

    //
    // 7. Run Intersect Bed
    //

    bedfilterbam( best, bedfiles )
    best = bedfilterbam.out.bam
    ch_versions = ch_versions.mix(bedfilterbam.out.versions)

    //
    // 8. Run Deamination workflow
    //

    deamination_stats( best, deduped.fixed )

    // get the meta-table from the "best"-libraries
    best = deamination_stats.out.best.map{ it[0] }
    fixed = deamination_stats.out.fixed.map{ it[0] }

    ch_final.mix( best ).mix( fixed ).set{ch_final}

    //
    // 9. Write the output files
    //

    write_reports( ch_final, ch_versions )
}
