include { SAMTOOLS_FILTER }  from '../modules/local/samtools_filter'
include { BAM_LENGTHFILTER } from '../modules/local/bam_lengthfilter'
include { SAMTOOLS_COUNT }   from '../modules/local/samtools_count'

workflow bamfilter {
    take:
        bam

    main:
        //
        // filter the bam files by samtools flag
        //

        SAMTOOLS_FILTER( bam )
        versions = SAMTOOLS_FILTER.out.versions.first()

        filtered = SAMTOOLS_FILTER.out.bam.map { meta, bam, count ->
            [ meta+[
                'ReadsRaw':count.text.split(',')[0].trim() as int,
                'ReadsFiltered':count.text.split(',')[1].trim() as int
            ], bam
            ]
        }

        //
        // filter the bam files by length
        //

        BAM_LENGTHFILTER( filtered )
        versions = versions.mix(BAM_LENGTHFILTER.out.versions.first())

        SAMTOOLS_COUNT( BAM_LENGTHFILTER.out.bam )
        versions = versions.mix( SAMTOOLS_COUNT.out.versions.first() )

        counted = SAMTOOLS_COUNT.out.bam.map { meta, bam, count ->
            [ meta+['ReadsLengthfiltered': count as int], bam ]
        }
    emit:
        bam = counted
        versions = versions
}