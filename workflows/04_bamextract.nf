include { GATHER_TAXON }  from "../modules/local/gather_taxon"
include { EXTRACT_TAXON }   from "../modules/local/extract_bam"
include { SAMTOOLS_SORT } from "../modules/local/samtools_sort"

workflow bamextract {
    take: bam
    take: translate

    main:

        //
        // 1. Create a list of reads for each taxon that need to be extracted
        //

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
        // Note: It will extract _more_ taxa than mapped, because the mapped taxa are stricter filtered than the extracted taxa
        // We dont have kmer-information here to filter on, only nuber of reads...
        EXTRACT_TAXON( ids )
        versions = EXTRACT_TAXON.out.versions.first()

        //
        // 3. Sort the bam
        //

        SAMTOOLS_SORT( EXTRACT_TAXON.out.bam )
        versions = versions.mix(SAMTOOLS_SORT.out.versions.first())

    emit:
        bam = SAMTOOLS_SORT.out.bam
        versions = versions

}

