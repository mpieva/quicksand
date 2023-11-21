// include modules that are used by the main workflow
include { MAP_BWA         } from '../modules/local/bwa'
include { SAMTOOLS_FILTER } from '../modules/local/samtools_filter'
include { SAMTOOLS_SORT   } from '../modules/local/samtools_sort'

standard_run = params.rerun ? false : true

workflow mapbam {
    take: bwa_in
    take: genomesdir
    take: ch_fixed
    take: ch_rerun

    main:
        // check if rerun,
        if( standard_run ){

            // prepare ch_fixed
            ch_fixed.map{ info ->
                [ info.Taxon, info ]
            }
            .set{ ch_fixed }

            //create a list of families that need replacement
            def famList = []
            specmap = new File("${params.fixed}")
            if(specmap.exists()){
                specmap.eachLine{famList << it.split("\t").flatten()[0]}
            }

            //now split the main channel into best and fixed
            bwa_in
            .branch{
                fixed: it[0].Family in famList
                best: true
            }
            .set{bwa_in}

            // prepare main channel
            // replace assignments in the fixed branch
            bwa_in.fixed
            .map{ meta, bam ->
                [meta.Family, meta, bam ]
            }
            .combine( ch_fixed, by:0 )
            .map{ fam, meta, bam, info -> //this is where the replacement happens
                [
                    meta+['Species':info.Species, "Reference":"fixed"], bam, file(info.Genome)
                ]
            }
            .set{fixed}

            // For the best-branch
            // prepare the genomesdir to match family/species

            genomesdir.map{ genome ->
                def fam = genome.getParent().getName()
                def matcher = genome.name =~ /(.+)\.fasta/
                def speciesName = matcher ? matcher[0][1] : null
                [[speciesName, fam], genome]
            }.set{ genomesdir }

            //combine with bam
            //prepare key (species,family)
            bwa_in.best.map{ meta, bam ->
                [
                    [meta.Species, meta.Family],
                    meta,
                    bam
                ]
            }
            .combine(genomesdir, by:0)
            .map{ key, meta, bam, genome ->
                [meta, bam, genome]
            }
            .set{ best }

            best.mix(fixed).set{ combined }
        } else {
            combined = ch_rerun
        }

        // Map with BWA
        MAP_BWA( combined )
        versions = MAP_BWA.out.versions.first()

        // Filter for Mapping quality
        SAMTOOLS_FILTER ( MAP_BWA.out.bam )
        versions = versions.mix( SAMTOOLS_FILTER.out.versions.first() )

        filtered = SAMTOOLS_FILTER.out.bam.map {
            [
              it[0]+['ReadsMapped':it[2].text.split(',')[1].trim() as int],
              it[1]
            ]
        }

        // Sort the bams
        SAMTOOLS_SORT( filtered )
        sorted = SAMTOOLS_SORT.out.bam
        versions = versions.mix( SAMTOOLS_SORT.out.versions.first() )

    emit:
        bam      = sorted
        versions = versions

}