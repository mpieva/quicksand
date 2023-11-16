include { BAM_RMDUP } from '../modules/local/bam_rmdup'


workflow dedupbam {
    take: mapped_bam
    main:
        // Take the mapped bam files
        // remove deduplicates

        BAM_RMDUP( mapped_bam )

        // Add the Number of Unique Sequences to the meta
        // first parse the bam-rmdup output file
        stats = BAM_RMDUP.out.txt.splitCsv(sep:'\t', header:true, limit:1)

        //now add the stats to the meta
        BAM_RMDUP.out.bam
        .combine(stats, by: 0)
        .map{ meta, bam, stats ->
            [
                meta+["ReadsDeduped":stats["out"] as int],
                bam
            ]
        }
        .set{ bam }

        versions = BAM_RMDUP.out.versions.first()



    emit:
        bam = bam
}