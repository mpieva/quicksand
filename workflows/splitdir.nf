include { ESTIMATE_CC }     from '../modules/local/ccestimate'
include { SAMTOOLS_FQ2BAM } from '../modules/local/samtools_fq2bam'

// some required functions
def has_ending(file, extension){
    return extension.any{ file.toString().toLowerCase().endsWith(it) }
}

workflow splitdir {
    take:
        split
    main:
        split
        .map{ [[:], it] }
            .branch {
                bam: it[1].getExtension() == 'bam'
                fastq: has_ending( it[1], ["fastq","fastq.gz","fq","fq.gz"])
                stats:  it[1].name =~ /split.*stats/
                fail: true
            }
            .set{ split }

        // Estimate cross-contamination if file exists
        ESTIMATE_CC( split.stats.first() )

        // convert fastq to bam
        SAMTOOLS_FQ2BAM( split.fastq )


    emit:
        bams = split.bam.mix( SAMTOOLS_FQ2BAM.out.bam )
        cc = ESTIMATE_CC.out.txt
        versions = SAMTOOLS_FQ2BAM.out.versions.first()
}