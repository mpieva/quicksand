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

        ch_bedfiles.map{ bedfile ->
            def matcher = bedfile.name =~ /(.+)\.masked\.bed/
            def speciesName = matcher ? matcher[0][1] : null
            [speciesName, bedfile]
        }
        .set{ ch_bedfiles }

        // now combine with the bedfile
        // so that bedtools intersect works
        // first, create the species as key
        ch_best.map{ meta, bam ->
            [ meta.Species, meta, bam ]
        }
        .combine( ch_bedfiles, by:0 ) // [sp, meta, bam, bed] --> now remove the sp again
        .map{ sp, meta, bam, bed ->
            [ meta, bam, bed ]
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
        .set{ ch_best }

    emit:
        bam = ch_best
        versions = versions
}