include { BEDTOOLS_INTERSECT } from '../modules/local/bedtools_intersect'
include { SAMTOOLS_COUNT     } from '../modules/local/samtools_count'
include { SAMTOOLS_COVERAGE  } from '../modules/local/samtools_coverage'
include { RUN_DUSTMASKER     } from '../modules/local/run_dustmasker'  
include { WRITE_BEDFILES     } from '../modules/local/write_bedfile'

workflow bedfilterbam {
    take: ch_best
    take: ch_bedfiles
    take: ch_fixed_references
    take: ch_fixed_mapped
    main:
        
        ch_new_bedfiles = Channel.empty()

        if (params.fixed_bedfiltering) {
            // with fixed references, we dont have the bedfile in the "masked" directory, so we need to run 
            // dustmasker (as we do in quicksand-build) and write the bedfile for each fixed reference

            ch_fixed_references.map{ meta ->
                [
                    meta, file(meta.Genome)
                ]
            }.set{ch_fixed_references}

            RUN_DUSTMASKER(ch_fixed_references)

            ch_acclist = RUN_DUSTMASKER.out.txt

            WRITE_BEDFILES(ch_acclist)

            ch_fixed_references = WRITE_BEDFILES.out.bed

            // now prepare the channel for later
            // mixing with ch_bedfiles
            ch_new_bedfiles = ch_fixed_references.map{ meta, bed ->
                    [meta.Species, bed]
                }
                .combine(          
                    ch_fixed_mapped.map{ meta, bam ->
                        [meta.Species, meta, bam]
                    }, by: 0 
                )
                .map{ sp, bed, meta, bam ->
                    [meta, bam, bed]
                }
        }
        
        // Take the "best subset"
        // find the species bed file in the masked folder

        // collect all bedfiles into one map [species: bedfile]
        // so that each entry can check whether its bedfile exists
        ch_bedfiles.map{ bedfile ->
            def matcher = bedfile.name =~ /(.+)\.masked\.bed/
            def speciesName = matcher ? matcher[0][1] : null
            [speciesName, bedfile]
        }
        .filter{ sp, bed -> sp }
        .toList()
        .map{ list -> list.collectEntries() }
        .set{ ch_bedfiles }

        // split into entries with and without a bedfile
        // entries without a bedfile skip the bedfiltering (instead of being dropped)
        ch_best.combine( ch_bedfiles )
        .branch{ meta, bam, beds ->
            bed: beds.containsKey(meta.Species)
            nobed: true
        }
        .set{ ch_best }

        ch_nobed = ch_best.nobed.map{ meta, bam, beds ->
            log.warn "[quicksand]: No bedfile '${meta.Species}.masked.bed' found for ${meta.id}. Skip bedfiltering"
            [ meta, bam ]
        }

        ch_best.bed.map{ meta, bam, beds ->
            [ meta, bam, beds[meta.Species] ]
        }
        .mix(ch_new_bedfiles) // its empty or the bedfiles of the fixed references
        .set{ ch_best }

        // And run bedtools instersect
        BEDTOOLS_INTERSECT( ch_best )

        ch_best = BEDTOOLS_INTERSECT.out.bam
        versions = BEDTOOLS_INTERSECT.out.versions.first()

        // get the counts
        SAMTOOLS_COUNT( ch_best )
        // add the counts to the meta
        ch_best = SAMTOOLS_COUNT.out.bam.map { meta, bam, count ->
            [ meta+[ 'ReadsBedfiltered': count as int ], bam ]
        }

        // get the coverage
        SAMTOOLS_COVERAGE( ch_best )
        // And add the coverage to the meta
        SAMTOOLS_COVERAGE.out.bam
        .map{ meta, bam, cov ->
            def (covered_bases, breadth, coverage) = cov.split()
            [
                meta+[ "PostBedCoveredBP":covered_bases.trim() as int ],
                bam
            ]
        }
        .mix( ch_nobed ) // ReadsBedfiltered and PostBedCoveredBP stay '-'
        .set{ ch_best }

    emit:
        bam = ch_best
        versions = versions
}