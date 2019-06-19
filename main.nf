#!/usr/bin/env nextflow

params.bam = ''
params.rg = ''
params.db = '/mnt/ramdisk/refseqReleaseKraken'
params.cutoff = 35
params.dedup = false    // deduplicate after splitting
params.level = 0        // bgzf compression level for intermediate files, 0..9


process splitBam {
    conda "$baseDir/envs/sediment.yaml"

    input:
    file 'input.bam' from file(params.bam)
    file 'indices.tsv' from file(params.rg)

    output:
    file 'split/*.bam' into splitfiles, split_to_extract mode flatten

    script:
    """
    mkdir split
    splitbam -c $params.level -d split -f indices.tsv --minscore 10 --maxnumber 0 input.bam
    """
}

splitfiles
    .map { [it.baseName, it] }
    .set { splitfiles }

split_to_extract
    .map { [it.baseName, it] }
    .set { split_to_extract }

(dedup_input, tofasta_input) = ( params.dedup
                             ? [splitfiles, Channel.empty()]
                             : [Channel.empty(), splitfiles] )

process removeDups {
    conda "$baseDir/envs/sediment.yaml"

    input:
    set rg, 'input.bam' from dedup_input
    
    output:
    set rg, 'dedup.bam' into dedup_to_assign, dedup_to_extract

    script:
    // rg = "${input_bam.baseName}"
    """
    countdups.py -l $params.level -o dedup.bam -s stat.txt -c $params.cutoff input.bam
    """
}

process toFasta {
    conda "$baseDir/envs/bam2fasta.yaml"

    input:
    set rg, 'input.bam' from tofasta_input.mix(dedup_to_assign)

    output:
    set rg, 'output.fa' into dedup_fasta

    script:
    """
    bam2fastx -a -Q -A -o output.fa input.bam
    """
}

process runKraken {
    conda "$baseDir/envs/sediment.yaml"
    publishDir 'kraken', mode: 'link'

    input:
    set rg, 'input.fa' from dedup_fasta

    output:
    set rg, "${kraken_translate}" into kraken_assignments
    set rg, "${kraken_out}" into kraken_raw

    script:
    kraken_out = "${rg}.kraken"
    kraken_translate = "${rg}.translate"
    """
    kraken -db $params.db --output $kraken_out input.fa
    kraken-translate -db $params.db --mpa-format $kraken_out >$kraken_translate
    """
}

process statsKraken {
    conda "$baseDir/envs/sediment.yaml"
    publishDir 'kraken', mode: 'link', saveAs: { "${rg}.report" }


    input:
    set rg, 'input.kraken' from kraken_raw

    output:
    file 'kraken.report'

    script:
    """
    kraken-mpa-report --db $params.db input.kraken >kraken.report
    """
}

for_extraction = ( params.dedup ? dedup_to_extract : split_to_extract )

// - Match split/deduped bam to kraken output based on readgroup
// - Filter for families under Mammalia
// - Emit unique (rg, bamfile, kraken_translate_file, family_name)
//
for_extraction
    .cross(kraken_assignments)
    .map { [it[0][0], it[0][1], it[1][1], it[1][1].readLines()] }
    .transpose()
    .filter { it[3] =~ /c__Mammalia.*f__/ }
    .map { rg, bam, kraken, asn -> [rg, bam, kraken, (asn =~ /f__([^|]*)/)[0][1]] }
    .unique()
    .set { for_extraction }

process extractBam {
    conda "$baseDir/envs/sediment.yaml"
    publishDir 'out/blast', mode: 'link', saveAs: { out_bam }

    input:
    set rg, 'input.bam', 'kraken.translate', family from for_extraction

    output:
    file 'output.bam' into extracted_reads
    set rg, family, stdout into extracted_read_count

    script:
    out_bam = "${family}/${rg}_extractedReads-${family}.bam"
    """
    extract_bam.py -f $family -k kraken.translate -c $params.level -o output.bam input.bam
    """
}

extracted_read_count
    .collectFile(storeDir: 'out/blast') { rg, family, count ->
        [ "${rg}.tsv", "${family}\t${count}"]
    }
