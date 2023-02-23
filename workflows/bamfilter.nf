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

        filtered = SAMTOOLS_FILTER.out.bam.map {
            [ it[0]+[
                'ReadsRaw':it[2].text.split(',')[0].trim() as int,
                'ReadsFiltered':it[2].text.split(',')[1].trim() as int
            ], it[1]
            ]
        }

        //
        // filter the bam files No. 2
        //

        BAM_LENGTHFILTER( filtered )
        versions = versions.mix(BAM_LENGTHFILTER.out.versions.first())

        SAMTOOLS_COUNT( BAM_LENGTHFILTER.out.bam )
        versions = versions.mix( SAMTOOLS_COUNT.out.versions.first() )

        counted = SAMTOOLS_COUNT.out.bam.map {
            [ it[0]+['ReadsLengthfiltered': it[2] as int], it[1] ]
        }
    emit:
        bam = counted
        versions = versions
}