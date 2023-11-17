// include modules that are used by the main workflow
include { MAP_BWA         } from '../modules/local/bwa'
include { SAMTOOLS_FILTER } from '../modules/local/samtools_filter'
include { SAMTOOLS_SORT   } from '../modules/local/samtools_sort'

workflow mapbam {
    take: bwa_in
    take: genomesdir
    main:
        // first, find the reference
        // prepare the genomesdir to match family/species

        genomesdir.map{ genome ->
            def fam = genome.getParent().getName()
            def matcher = genome.name =~ /(.+)\.fasta/
            def speciesName = matcher ? matcher[0][1] : null
            [[speciesName, fam], genome]
        }.set{ genomesdir }

        //combine with bam
        //prepare key (species,family)
        bwa_in.map{ meta, bam ->
            [
                [meta.Species, meta.Family],
                meta,
                bam]
        }
        .combine(genomesdir, by:0)
        .map{ key, meta, bam, genome ->
            [meta, bam, genome]
        }
        .set{ bwa_in }

        // Map with BWA
        MAP_BWA( bwa_in )
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