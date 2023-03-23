include { GATHER_TAXON }  from "../modules/local/gather_taxon"
include { EXTRACT_BAM }   from "../modules/local/extract_bam"
include { SAMTOOLS_SORT } from "../modules/local/samtools_sort"

workflow bamextract {
    take: bam
    take: translate

    main:
        versions = channel.empty()

        //
        // 1. Create a list of reads for each taxon that need to be extracted
        //

        // parse the translate-file
        translate.map{ meta, translate ->
            [meta, translate, translate.readLines()]
        }
        .transpose()
        .filter{ it[2] =~ "${params.taxlvl}__" }
        .map{ meta, translate, asgn ->
            [meta + ['taxon':(asgn =~ "${params.taxlvl}__([^|]*)")[0][1]],
            translate]
        }
        .unique()
        .set{translate}

        GATHER_TAXON(translate)

        //
        // 2. Use the list of reads to extract a subset of bams of the respective taxon
        //

        // prepare the bam file
        bam.map{ meta, bam ->
            [meta.id, bam]
        }
        .set{ bam }

        // remodel gather-output channel and combine with sorted bams
        GATHER_TAXON.out.ids
        .map{ meta, ids, count ->
            [ meta.id, meta + ['ReadsExtracted':count.trim() as int], ids ]
        }
        .combine(bam, by:0)
        .map{ id, meta, ids, bam ->
            [meta, ids, bam]
        }
        .set{ ids }

        // extract the bams
        EXTRACT_BAM( ids )
        versions = versions.mix(EXTRACT_BAM.out.versions.first())

        //
        // 3. Sort the bam
        //

        SAMTOOLS_SORT( EXTRACT_BAM.out.bam )
        versions = versions.mix(SAMTOOLS_SORT.out.versions.first())

    emit:
        bam = SAMTOOLS_SORT.out.bam
        versions = versions

}

