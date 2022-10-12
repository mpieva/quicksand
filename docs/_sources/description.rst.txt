Description
===========

.. image:: images/pipeline_overview_v1.6.png
	:width: 800
	:align: center
	:alt: Graphical overview over the processes of the quicksand pipeline


This section describes the processes implemented in the **quicksand** pipeline. See :ref:`quicksand_build-page` for a description of
the quicksand-build pipeline

Workflow
--------

The basic workflow of the quicksand pipline consists of a **metagenomic classification**, a **mapping** and a basic **analysis** of the alignment. The pipeline
produces taxonomic profiles for the analyzed input files, as well as the binned :file:`bam` files, grouped by redgroup and taxonomy for all levels of the
pipeline - as shown in the figure above. Each process is in detail described in the following chapters

EstimateCrossContamination
""""""""""""""""""""""""""
WIP

fastq2bam
""""""""""

All :file:`fastq` files entering the pipeline are converted to unpaired and unmapped :file:`bam` files using :code:`samtools import -0`.
In case of paired end :file:`fastq` files, make sure to merge them before using the pipeline.

fiterBam
"""""""""

:file:`Bam` files are filtered based on the provided filterflag using :code:`samtools view -F`. By default paired-end reads are removed from the
:file:`bam` file.

filterLength
""""""""""""

Using a custom :code:`bam-lengthfilter` package, sequences are removed from the :file:`bam` file that fall below the :code:`--bamfilter_length_cutoff` threshold (default: 35bp)

toFasta
""""""""

As preparation for the krakenuniq run, :file:`bam` files are converted to :file:`fasta` files

runKrakenUniq
""""""""""""""

Sequences contained in the :file:`fasta` file are classified by the metagenomic classifier :code:`krakenuniq`. Quicksand uses krakenuniq with a precompiled
database created from the current non-redundant mtDNA RefSeq database with a default kmer size of 22 (see :ref:`quicksand_build-page` or :ref:`setup`).
The speed of krakenuniq allows for a quick sorting of sequences into families. To filter out false-positive assignments, families are
removed from the assignment falling below the minimum number of reads (:code:`--krakenuniq_min_reads`) and the minimal number of kmers (:code:`--krakenuniq_min_kmers`)
on the family node. The result of this process is a taxonomic profile of the readgroup.

findBestNode
""""""""""""

This process parses the kraken-reports. For each assigned family reported by krakenuniq, the node with the highest number of assigned unique kmers is
picked as the taxon representative for that family.

extractBam
""""""""""

Using the kraken-report and the length-filtered :file:`bam` file, this process collects all sequences assigned to one clade into a new :file:`bam`
file. Extraction happens either on the family or order-level, as specified with the :code:`--taxlvl` flag, using the custom :code:`bamfilter` package.

mapBWA
""""""

The extracted sequences are mapped against all the reference genomes of species belonging to the 'bestNode' found in the 'findBestNode' process
using the :code:`bwa bam2bam` command of the `network-aware fork <https://github.com/mpieva/network-aware-bwa>`_ of BWA with
ancient parameters (:code:`n 0.01 -o 2 -l 16500`). Unmapped sequences or sequences with a mapping quality of less than 25 are removed from the alignment


filterMappedBam
""""""""""""""""

Mapped bam files are filtered for the set alignment quality score

dedupBam
""""""""

Exact PCR-duplicates are collapsed into unique sequences using `bam-rmdup <https://github.com/mpieva/biohazard-tools>`_ based on the sharing of identical
alignment start and end coordinates. From all mapped genomes, the one with the highest numbers of basepairs covered is picked as _the_
representative species for the subsequent steps.

runIntersectBed
""""""""""""""""

The deduped alignments are then depleted of reads that overlap sites marked as non-informative by :code:`dustmasker`. That step is skipped
for families with a fixed reference genome (see :code:`--fixed` flag)

analyzeDeamination
""""""""""""""""""

The final step(s) in the pipeline look for C to T substitutions in the query sequences in respect to the aligned reference genome.
Ancient DNA shows characteristic C to T substitutions at the 3’ and 5’ ends
of DNA fragments - a degradation pattern used to identify ancient DNA. Families which sequences show more than 10% of terminal C bases in the
reference genome replaced by a T are reported as being ancient (++).

extractDeaminatedReads
""""""""""""""""""""""

For families with a fixed reference genome (see :code:`--fixed` flag), extract deaminated reads into two different :code:`.bam` files. One file with sequences that
show C to T substitutions on one of the terminal base pairs and one file that looks at the first and last 3 terminal base positions.

maskDeamination
""""""""""""""""

For families with a fixed reference genome (see :code:`--fixed` flag), take the extracted deaminated reads and set the alignment quality score of the last
three terminal T bases to 15.

createMpileups
""""""""""""""

