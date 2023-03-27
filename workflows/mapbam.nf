// include modules that are used by the main workflow
include { MAP_BWA         } from '../modules/local/bwa'
include { SAMTOOLS_FILTER } from '../modules/local/samtools_filter'

workflow mapbam {
    take: bwa_in
    main:

        MAP_BWA( bwa_in )
        versions = MAP_BWA.out.versions.first()

        // filter mapped bams (quality)
        // sort mapped bams
        // count mapped bams



    emit:
        versions = versions

}