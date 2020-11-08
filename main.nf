#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    Usage: nextflow run sediment_nf --bam PATH --rg PATH --db PATH --genome PATH --bedfiles
                                    [--cutoff N] [--level N] [--keeppaired] [--filterunmapped] [--krakenfilter N]
           
    Run sediment analysis pipeline.
    
    required arguments:
      --bam PATH         input BAM file
      --rg  PATH         tab-separated file containing index combinations
                         - format of file: 'LibID<tab>P7<tab>P5'
      --db PATH          Kraken database (set a default in nextflow.config)
      --genome PATH      genome for alignment (set a default in nextflow.config)
      --bedfiles PATH    bed-files masking the genomes for alignment (set a default in nextflow.config)
      
    optional arguments:
      --cutoff N         length cutoff (default: 35)
      --quality N        quality filter (default: 25)
      --bwacutoff N      cutoff for number of kraken matches (default: 0)
      --keeppaired       keep paired reads (default: filter paired)
      --filterunmapped   filter unmapped reads
      --krakenfilter N   kraken-filter with threshold N [0,1] (default: 0)
      --level N          set BGZF compression level (default: 6)
      --krakenthreads N  numbrer of threads per Kraken process (default: 4)
      
    A selection of built-in Nextflow flags that may be of use:
      -resume            resume processing; do not re-run completed processes
      -profile NAME      use named execution profile
      -qs N              queue size; number of CPU cores to use (default: all)
      -N EMAIL           send completion notifcation to email address
    
    The following execution profiles are available (see '-profile' above):
      * standard         execute all processes on local host
      * cluster          execute certain CPU-intenisve processes on SGE cluster
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

//All params can be found in the config-file

process splitBam {
    conda "$baseDir/envs/sediment.yaml"
    maxForks 1
    publishDir 'split', mode: 'link'
    label 'local'

    input:
    file 'input.bam' from file(params.bam)
    file 'indices.tsv' from file(params.rg)

    output:
    file '*.bam' into splitfiles mode flatten
    stdout into splitscriptstats

    script:
    """
    splitbam -s -c $params.level -f indices.tsv --minscore 10 --maxnumber 0 input.bam
    """
}


splitscriptstats
    .collectFile(storeDir: 'stats', name: "splitstats.tsv", newLine: true)

splitfiles
    .map{ [it.baseName, it] }      //it.basename = rg, it=full bam file
    .into{ splitfiles; splitstats }


process splitStats {
    conda "$baseDir/envs/sediment.yaml"
    tag "$rg"

    input:
    set rg,'input.bam' from splitstats

    output:
    set rg, stdout into splitcounts

    script:
    """
    samtools view -c input.bam
    """
}

// if keeppaired==True, use an empty channel, else use splitfiles

filter_paired_in = params.keeppaired ? Channel.empty() : splitfiles

process filterPaired {
    conda "$baseDir/envs/sediment.yaml"
    tag "$rg"

    input:
    set rg, 'input.bam' from filter_paired_in

    output:
    set rg, 'output.bam' into filter_paired_out

    script:
    """
    samtools view -b -u -F 1 -o output.bam input.bam
    """
}
//here the paths come together again
post_filter_paired = params.keeppaired ? splitfiles : filter_paired_out

// and do the same with the filter-unmapped step
filter_unmapped_in = params.filterunmapped ? post_filter_paired : Channel.empty()

process filterUnmapped {
    conda "$baseDir/envs/sediment.yaml"
    tag "$rg"

    input:
    set rg, 'input.bam' from filter_unmapped_in

    output:
    set rg, 'output.bam' into filter_unmapped_out

    script:
    """
    samtools view -b -u -F 4 -o output.bam input.bam
    """
}

post_filter_unmapped = params.filterunmapped ? filter_unmapped_out : post_filter_paired

process filterLength {
    conda "$baseDir/envs/sediment.yaml"
    tag "$rg"

    input:
    set rg, 'input.bam' from post_filter_unmapped

    output:
    set rg, 'output.bam', stdout into tofasta_in
    set rg, 'output.bam' into for_extraction
    set rg, stdout into filtercounts

    script:
    """
    bam-lengthfilter -c $params.cutoff -l $params.level -o output.bam input.bam
    samtools view -c output.bam
    """
}

filtercounts.join(splitcounts)
    .map { rg, fc, sc -> "${rg}\t${sc.trim()}\t${fc.trim()}"}
    .collectFile(storeDir: 'stats', name: "splitcounts.tsv", newLine: true,
                 seed: "readgroup\tsplit count\tfiltered count")

process toFasta {
    conda "$baseDir/envs/sediment.yaml"
    tag "$rg"

    input:
    set rg, 'input.bam', filtered_count from tofasta_in

    output:
    set rg, 'output.fa' into tofasta_out

    when:
    filtered_count.toInteger() > 0

    script:
    """
    samtools fasta input.bam > output.fa
    """
}

process runKraken {
    conda "$baseDir/envs/sediment.yaml"
    publishDir 'kraken', mode: 'link', saveAs: {"${rg}.kraken"}
    cpus "${params.krakenthreads}"
    memory '16GB'
    label 'bigmem'
    label 'local'
    tag "$rg"

    input:
    set rg, "input.fa" from tofasta_out

    output:
    set rg, "output.kraken" into kraken_out

    script:
    """
    kraken --threads ${task.cpus} --db $params.db --output output.kraken --fasta-input input.fa
    """
}

kraken_filter_in = params.krakenfilter ? kraken_out : Channel.empty()

process filterKraken{
    conda "$baseDir/envs/sediment.yaml"
    publishDir 'kraken', mode: 'link', saveAs: {"${rg}.kraken_filter"}
    label 'local'
    tag "$rg"

    input:
    set rg, "input.kraken" from kraken_filter_in

    output:
    set rg, "output_filtered.kraken" into kraken_filter_out

    when:
    params.krakenfilter

    script:
    """
    kraken-filter -threshold $params.krakenfilter --db $params.db input.kraken > output_filtered.kraken
    """
}

post_kraken_filter = params.krakenfilter ? kraken_filter_out : kraken_out

process translateKraken{
    conda "$baseDir/envs/sediment.yaml"
    publishDir 'kraken', mode: 'link', pattern:"*translate", saveAs: {"${rg}.translate"}
    publishDir 'kraken', mode: 'link', pattern:"*report", saveAs: {"${rg}.report"}
    label 'local'
    tag "$rg"

    input:
    set rg, "input.kraken" from post_kraken_filter

    output:
    set rg, "kraken.translate" into kraken_assignments
    set rg, "kraken.report" into find_best

    script:
    """
    kraken-translate --db $params.db --mpa-format input.kraken > kraken.translate
    kraken-report --db $params.db input.kraken > kraken.report
    """
}
process findBestSpecies{
    conda "$baseDir/envs/sediment.yaml"
    tag "$rg"

    input:
    set rg, "kraken.report" from find_best

    output:
    set rg, "parsed_record.tsv" into best_species

    script:
    """
    parse_report.py kraken.report
    """
}

Channel.fromPath("${params.genome}/taxid_map.tsv", type:'file')
    .splitCsv()
    .map{it[0].split('\t').flatten()}
    .map{[it[0], it[2]]}
    .groupTuple()
    .set{taxid}

best_species
    .map{[it[0], it[1].readLines()]}
    .transpose()
    .map{[it[1].split('\t')[1], it[0], it[1].split('\t')[0]]}
    .combine(taxid, by:0)  //[taxid, readgroup, family, [species, species]]
    .transpose()
    .map{[it[1]+it[2], it[3]]}
    .set{best_species}

for_extraction
    .cross(kraken_assignments)
    .map { [it[0][0], it[0][1], it[1][1], it[1][1].readLines()] }
    .transpose()
    .filter { it[3] =~ /c__Mammalia.*f__./ }
    .map { rg, bam, kraken, asn -> [rg, bam, kraken, (asn =~ /f__([^|]*)/)[0][1]] }
    .unique()
    .set { for_extraction }

process gatherByFamily {
    tag "$rg:$family"

    input:
    set rg, 'input.bam', 'kraken.translate', family from for_extraction

    output:
    set family, rg, 'input.bam', 'ids.txt', stdout into (prepared_for_extraction, count_for_stats)

    script:
    """
    grep "c__Mammalia.*f__$family" kraken.translate | cut -f1 | tee ids.txt | wc -l
    """
}

count_for_stats
        .collectFile(storeDir: 'stats') { family, rg, bamf, idf, count ->
            [ "${rg}_extracted.tsv", "${family}\t${count}"]
        }

process extractBam {
    conda "$baseDir/envs/sediment.yaml"
    publishDir 'out', mode: 'link', saveAs: {"${family}/${rg}_extractedReads-${family}.bam"}
    tag "$rg:$family"

    input:
    set family, rg, 'input.bam', 'ids.txt', idcount from prepared_for_extraction

    output:
    set rg, family, 'output.bam' into extracted_reads

    when:
    idcount.toInteger() >= params.bwacutoff

    script:
    """
    bamfilter -i ids.txt -l $params.level -o output.bam input.bam
    """
}

extracted_reads
    .map{[it[0]+it[1], it[0], it[1], it[2]]}
    //extracted reads --> [new_Key, rg, fam, ExtractedReads_Hominidae.bam]
    //best_species --> [New_key, species]
    .cross(best_species)
    .map{x,y -> [x[1], x[2], y[1], x[3]] }
    //[Readgroup, Hominidae, Homo_sapiens, ExtractedReads_Hominidae.bam]
    .set{extracted_reads}

process mapBwa {
    publishDir 'out', mode: 'link', saveAs: { out_bam }, pattern: '*.bam'
    conda "$baseDir/envs/sediment.yaml"
    tag "$rg:$family:$species"

    input:
    set rg, family, species, "input.bam" from extracted_reads

    output:
    set family, rg, species, 'output.bam' into mapped_bam
    set family, rg, species, stdout into mapped_count
    set family, rg, species, 'coverage.txt' into coverage_count

    script:
    out_bam = "${family}/aligned/${rg}.${species}.bam"
    """
    samtools sort -n -l0 input.bam \
    | $params.bwa bam2bam -g ${params.genome}/$family/\"${species}.fasta\"  -n 0.01 -o 2 -l 16500 --only-aligned - \
    | samtools view -b -u -q $params.quality \
    | samtools sort -l $params.level -o output.bam
    samtools coverage -H output.bam | cut -f 5 > coverage.txt
    samtools view -c output.bam
    """
}

mapped_count
        .collectFile(storeDir: 'stats') { family, rg, species, count ->
            [ "${rg}_mapped.tsv", "${family}\t${species}\t${count}"]
        }

coverage_count
        .collectFile(storeDir: 'stats') { family, rg, species, cov_file ->
            [ "${rg}_mapped_coverage.tsv", "${family}\t${species}\t" + cov_file.text]
        }

process dedupBam {
    publishDir 'out', mode: 'link', saveAs: {"${family}/aligned/${rg}.${species}_deduped.bam"}
    conda "$baseDir/envs/sediment.yaml"
    tag "$rg:$family:$species"

    input:
    set family, rg, species, 'input.bam' from mapped_bam

    output:
    set family, rg, species, stdout into deduped_count
    set species, family, rg, "output.bam", stdout into deduped_bam
    file 'output.bam'

    script:
    """
    $params.bamrmdup -r -o output.bam input.bam > rmdup.txt
    samtools view -c output.bam
    """
}
deduped_bam
	/*sort by readgroup, family and if same (?:) by count. toList() --> wait for all dedup-processes to finish
	flatMap --> Emit every record in the list separatly (for groupTuple)
	groupTuple by readgroup and family
	Map only the best species per family*/
	.toSortedList({
		a,b -> a[2]+a[1] <=> b[2]+b[1] ?: a[-1] as int <=> b[-1] as int
		})
	.flatMap{n -> n[0..-1]}
	.groupTuple(by:[2,1])	//[[sp,sp,sp], family, rg, [bam,bam,bam],[count < count < count]]
	.map{n -> [n[0][-1], n[1], n[2], n[3][-1]]} //[species, family, rg, bamfile]
	.set{best_deduped}


deduped_count
    .collectFile(storeDir: 'stats') { family, rg, species, count ->
    [ "${rg}_unique_mapped.tsv", "${family}\t${species}\t${count}"]
    }


Channel.fromPath("${params.bedfiles}/*.bed", type:'file')	//all the bedfiles
    .map{[it.baseName.replaceAll(/.masked/,""), it] }		//make a map, with species as key
    .cross(best_deduped)									//throw it together
    .map{x, y -> [y[0],y[1],y[2],y[3],x[1]]}				//get species, family, rg, bamfile and bedfile
    .set{to_bed}


//and filter out reads that intersect with masked regions
process runIntersectBed{
    tag "$rg:family:species"
    publishDir 'out', mode: 'link', saveAs: {"${family}/bed/${rg}.${species}_deduped_bedfiltered.bam"}
    conda "$baseDir/envs/sediment.yaml"

    input:
    set species, family, rg, "inbam.bam", "inbed.bed" from to_bed

    output:
    file "outbam.bam"
    set family, rg, species, stdout into bedfilter_count
        
    script:
    """
    bedtools intersect -a inbam.bam -b inbed.bed -v > outbam.bam
    samtools view -c outbam.bam
    """
    }
    
bedfilter_count
        .collectFile(storeDir: 'stats') { family, rg, species, count ->
            [ "${rg}_bedfiltered.tsv", "${family}\t${species}\t${count}"]
        }
