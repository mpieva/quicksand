.. _output-page:

Input and Output
================

.. _input:

Input
-----

quicksand requires as input demultiplexed, merged, and adapter trimmed sequencing libraries in :file:`.fastq` or :file:`.bam` format.
Use the :code:`--split` flag to point to the directory containing these files. quicksand refers to the name of the files
as readgroups::

    splitdir/
        readgroup1.fastq
        readgroup2.fastq
        readgroup3.bam

.. note::
    Quicksand will process all the :file:`.fastq,fq,fastq.gz,fq.gz` and :file:`.bam`-files within the :code:`--split` directory and ignore all other files.


Output
------

quicksand writes all output files to the **quicksand_v2.2** directory. Within this directory the files are
layed out as follows::

    quicksand_v2.2
    ├── out
    │    └── {taxon}
    │         ├── 1-extracted
    │         │    └── {RG}_extractedReads-{taxon}.bam
    │         ├── best // (for families not in --fixed)
    │         │    ├── 2-aligned
    │         │    │     └── {RG}.{family}.{species}.bam
    │         │    ├── 3-deduped
    │         │    │     └── {RG}.{family}.{species}_deduped.bam
    │         │    └── 4-bedfiltered
    │         │          └── {RG}.{family}.{species}_deduped_bedfiltered.bam
    │         └── fixed // (for families in --fixed)
    │              ├── 2-aligned
    │              │     └── {RG}.{family}.{species}.bam
    │              ├── 3-deduped
    │              │     └── {RG}.{family}.{species}_deduped.bam
    │              ├── 5-deaminated
    │              │     ├── {RG}.{family}.{species}_deduped_deaminated_1term.bam
    │              │     └── {RG}.{family}.{species}_deduped_deaminated_3term.bam
    │              └── 6-mpileups
    │                    ├── {RG}.{family}.{species}_term1_mpiled.tsv
    │                    ├── {RG}.{family}.{species}_term3_mpiled.tsv
    │                    └── {RG}.{family}.{species}_all_mpiled.tsv
    ├── stats
    │    ├── splitcounts.tsv
    │    ├── {RG}.kraken.report
    │    ├── {RG}.kraken.translate
    │    ├── {RG}_00_extracted.tsv
    │    ├── {RG}_01_mapped.tsv
    │    ├── {RG}_02_deduped.tsv
    │    ├── {RG}_03_bedfiltered.tsv
    │    └── {RG}_04_deamination.tsv
    ├── nextflow
    │    ├── {DATE}.commands
    │    └── {DATE}.config
    ├── work
    │    └── ...
    ├── cc_estimates.tsv
    ├── filtered_report_{N}p_{N}b.tsv
    └── final_report.tsv


.. _files:

Files
"""""

Directory: out/TAXON
~~~~~~~~~~~~~~~~~~~~

.. rst-class:: file
*1-extracted/$\{RG\}.extractedReads-$\{taxon\}.bam*

.. rst-class:: description
BAM FILE. Contains the DNA sequences of one readgroup assigned by KrakenUniq to one taxon [family or order].

.. rst-class:: file
*2-aligned/$\{RG\}.$\{family\}.$\{species\}.bam*

.. rst-class:: description
BAM FILE. Contains the aligend sequences after mapping the extractedReads to the reference species

.. rst-class:: file
*3-deduped/$\{RG\}.$\{family\}.$\{species\}_deduped.bam*

.. rst-class:: description
BAM FILE. The same alignment, but depleted of PCR duplicates.

.. rst-class:: file
*4-bedfiltered/$\{RG\}.$\{family\}.$\{species\}_deduped_bedfiltered.bam*

.. rst-class:: description
| BAM FILE. The deduped alignment, but depleted of reads overlapping low-complexity regions
| specified in the provided bedfiles for the given species.

.. rst-class:: file
*5-deaminated/$\{RG\}.$\{family\}.$\{species\}_deduped_deaminated_1term.bam*

.. rst-class:: description
| BAM FILE. The deduped alignment, filtered for reads that show a C to T
| substitution at one of the terminal positions in respect to the reference genome

.. rst-class:: file
*5-deaminated/$\{RG\}.$\{family\}.$\{species\}_deduped_deaminated_3term.bam*

.. rst-class:: description
| BAM FILE. The deduped alignment, filtered for reads that show a C to T
| substitution at one of the terminal `three` positions in respect to the reference genome

.. rst-class:: file
*6-mpileups/$\{RG\}.$\{family\}.$\{species\}_all_mpiled.tsv*

.. rst-class:: description
| TSV FILE. The deduped alignment, but in mpileup format.
| The first three positions of each sequence are masked by setting the mapping quality to 0

.. rst-class:: file
*6-mpileups/$\{RG\}.$\{family\}.$\{species\}_1term_mpiled.tsv*

.. rst-class:: description
| TSV FILE. Mpileup format. The first three positions of each sequence are masked by setting the mapping quality to 0.
| The pileup contains only reads showing a C to T substitution at one of the terminal positions in respect to the reference genome

.. rst-class:: file
*6-mpileups/$\{RG\}.$\{family\}.$\{species\}_3term_mpiled.tsv*

.. rst-class:: description
| TSV FILE. Mpileup format. The first three positions of each sequence are masked by setting the mapping quality to 0.
| The pileup contains only reads showing a C to T substitution at one of the terminal `three` positions in respect to the reference genome


Directory: stats
~~~~~~~~~~~~~~~~

.. rst-class:: file
*$\{RG\}.report*

.. rst-class:: description
The standard krakenuniq report

.. rst-class:: file
*$\{RG\}.translate*

.. rst-class:: description
The human readable kraken report in mpa-format

.. rst-class:: file
*stats/splitcounts.tsv*

.. rst-class:: description
TSV FILE. Contains for each readgroup the number of reads before (raw) and after the initial filter step::

    RG          ReadsRaw      ReadsFiltered ReadsLengthfiltered
    test1       235           235           230
    test2       235           235           230
    test3       235           235           230

.. rst-class:: file
*$\{RG\}_00_extracted.tsv*

.. rst-class:: description
TSV FILE. Contains the number of sequences assigned to a taxon based on the KrakenUniq classification::

    Taxon       ReadsExtracted
    Hominidae   235

.. rst-class:: file
*$\{RG\}_01_mapped.tsv*

.. rst-class:: description
TSV FILE. Contains for each readgroup and family the number of sequences mapped to the reference genome. The column 'Reference' shows if the reference
genome was fixed. The proportion mapped is the proportion of mapped to extracted reads::

    Order     Family      Species       Reference    ReadsMapped   ProportionMapped
    Primates  Hominidae   Homo_sapiens  fixed        235           0.913

.. rst-class:: file
*$\{RG\}_02_deduped.tsv*

.. rst-class:: description
TSV FILE. Contains for each readgroup and family the number of unique reads mapped to the reference genome, the duplication rate
and information from the :code:`samtools coverage` command::

    Order: The taxonomic order
    Family: The taxonomic family
    Species: The taxonomic species used as reference for mapping
    Reference: The reference type: either 'best' or 'fixed'
    ReadsDeduped: The number of unique reads
    DuplicationRate: The duplication rate of the unique reads
    CoveredBP: 'covbases' of the samtools coverage command: The number of covered bases in the reference genome
    Coverage: 'meandepth' of the samtools coverage command: The mean depth of coverage
    Breadth: 'coverage' of the samtools coverage command (by 100): the proportion of covered bases in the reference genome
    ExpectedBreadth: Expected breadth based on the inStrain formula: expected_breadth = 1-e^(-0.833*coverage). See
        https://instrain.readthedocs.io/en/latest/important_concepts.html
    ProportionExpectedBreadth: The proportion of Breadth / ExpectedBreadth

.. rst-class:: file
*$\{RG\}_03_bedfiltered.tsv*

.. rst-class:: description
TSV FILE. Contains for each readgroup and family the number of sequences remaining in the bam-file after bedfiltering and the number of covered basepairs
in the reference genome after removal of low-complexity sequences::

    Order     Family      Species       Reference  ReadsBedfiltered PostBedCoveredBP
    Primates  Hominidae   Homo_sapiens  fixed      97               4177

.. rst-class:: file
*$\{RG\}_04_deamination.tsv*

.. rst-class:: description
TSV FILE. Contains for each readgroup the deamination stats for the BAM file after bedfiltering::

          Ancientness:  ++  = more than 9.5% of the reads that show a terminal C in both the 5' and 3' position in the reference genome, carry a T
                        +   = more than 9.5% of the reads that show a terminal C in either the 5' or 3' position in the reference genome, carry a T
                        -   = no signs for DNA deamination patterns

          ReadsDeam(1term): The number of reads (after deduplication and bedfiltering) that show a deamination in the terminal base positions
          ReadsDeam(3term): The number of reads (after deduplication and bedfiltering) that show a deamination in the three terminal base positions
          Deam5(95ci):      For the terminal 5' end, the percentage of C to T substitutions (and the 95% confidence interval)
          Deam3(95ci):      For the terminal 3' end, the percentage of C to T substitutions (and the 95% confidence interval)
          Deam5Cond(95ci):  Taken only 3' deaminated sequences, report the percentage of C to T substitutions (and the 95% confidence interval) at the 5' terminal base
          Deam3Cond(95ic):  Taken only 5' deaminated sequences, report the percentage of C to T substitutions (and the 95% confidence interval) at the 3' terminal base


final_report.tsv
~~~~~~~~~~~~~~~~

The final report contains all the columns presented above. In Addition, the final report contains a column :code:`FamPercentage` which provides the relative
proportion of *final reads* (after deduplication or bedfiltering) of the assigned family in the readgroup. If there are several lines for one family and readgroup (e.g. after a rerun or multiple fixed references)
the highest number of final reads is used as the baseline for the other entries of the same family

filtered_report.tsv
~~~~~~~~~~~~~~~~

The filtered report contains all the columns from the final_report. However, the report is filtered by the two values :code:`FamPercentage` and :code:`ProportionExpectedBreadth` as
provided by the flags :code:`--reportfilter_percentage` and :code:`--reportfilter_breadth` (both default to 0.5).


| The :file:`cc_estimates.tsv` files contains information about index-hopping and cross contamintaion
| The :file:`nextflow` directory contains information about the run, like the commandline used and the config-files provided
| the :file:`work` directory can be deleted after the run - it contains nextflow specific intermediate files
