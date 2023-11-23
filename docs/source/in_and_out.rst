.. _output-page:

.. _output:

Input and Output
================

.. _input:

Input
-----

The pipeline uses as input :file:`.fastq` or :file:`.bam` files that contain demultiplexed, merged, and adapter trimmed reads.
Use the :code:`--split` flag to point to the directory that contains these files. The pipeline refers to the name of the files
as readgroups. The reads within the files are assigned, processed and structured by readgroups.
:file:`.bam` and :file:`.fastq` files can be mixed::

    splitdir/
        readgroup1.fastq
        readgroup2.fastq
        readgroup3.bam

.. note::
    Quicksand will process all the :file:`.fastq,fq,fastq.gz,fq.gz` and :file:`.bam`-files within the :code:`--split` directory and ignore all remaining files.
    Be sure to name/rename your files accordingly

.. _output:

Output
------

Several directories and files should appear after the run. The 'taxon' corresponds to either the family or the order level name::

    quicksand_v2.0
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
    └── final_report.tsv


.. _files:

Files explained
"""""""""""""""

The content of the files is explained here:

**out/**

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - File
     - Description
   * - :file:`$\\{RG\\}.extractedReads-$\\{taxon\\}.bam`
     - Contains all DNA sequences from one readgroup assigned by krakenuniq to one taxon (family or order).
   * - :file:`aligned/$\\{RG\\}.$\\{family\\}.$\\{species\\}.bam`
     - An alignemnt file. The result of mapping all extractedReads (see above) against the genome of the assigned 'best' species
   * - :file:`aligend/$\\{RG\\}.$\\{family\\}.$\\{species\\}_deduped.bam`
     - The same alignment file, but depleted of PCR duplicates - reads with the same start and end coordinates in the alignment
   * - :file:`bed/$\\{RG\\}.$\\{family\\}.$\\{species\\}_deduped_bedfiltered.bam`
     - The same aligned and deduped bamfile, but additionally depleted of reads overlapping the low-complexity regions specified in the :code:`--bedfiles` for the given species
   * - :file:`deaminated/$\\{RG\\}.$\\{family\\}.$\\{species\\}_deduped_deaminated_1term.bam`
     - | The aligned and deduped bamfile, filtered for reads that show a C to T
       | substitution at one of the terminal basepair positions in respect to the reference genome
   * - :file:`deaminated/$\\{RG\\}.$\\{family\\}.$\\{species\\}_deduped_deaminated_3term.bam`
     - | The aligned and deduped bamfile, filtered for reads that show a C to T
       | substitution at one of the terminal three basepair positions in respect to the reference genome
   * - :file:`deaminated/$\\{RG\\}.$\\{family\\}.$\\{species\\}_all_mpiled.tsv`
     - The aligned and deduped bamfile, but in mpileup format
   * - :file:`deaminated/$\\{RG\\}.$\\{family\\}.$\\{species\\}_1term_mpiled.tsv`
     - The deaminated 1term bamfile, with masked terminal T bases - in mpileup format
   * - :file:`deaminated/$\\{RG\\}.$\\{family\\}.$\\{species\\}_3term_mpiled.tsv`
     - The deaminated 3term bamfile, with masked terminal T bases - in mpileup format


**stats/**

.. list-table::
   :widths: 20 80
   :header-rows: 1

  * - File
    - Description
  * - :file:`$\\{RG\\}.report`
    - The standard krakenuniq report
  * - :file:`$\\{RG\\}.translate`
    - The read-wise human readable kraken report in mpa-format
  * - :file:`stats/splitcounts.tsv`
    - | For each readgroup (RG), show the number of reads before (raw) and after the filterBam process, as
      | well as the number of reads after the bam-lengthfilter process
      ::

          RG          ReadsRaw      ReadsFiltered ReadsLengthfiltered
          test1       235           235           230
          test2       235           235           230
          test3       235           235           230
  * - :file:`$\\{RG\\}_00_extracted.tsv`
    - Shows the number of reads extracted for each assigned taxon based on the kraken assignments
      ::

          Taxon       ReadsExtracted
          Hominidae   235
  * - :file:`$\\{RG\\}_01_mapped.tsv`
    - | For each readgroup (RG) show the number of reads mapped against the reference genome, if the reference genome was fixed
      | (see :code:`--fixed` flag) and the proportion of mapped reads (from the number of extracted reads for this family)
      ::

          Order     Family      Species       Reference    ReadsMapped   ProportionMapped
          Primates  Hominidae   Homo_sapiens  fixed        235           0.913
  * - :file:`$\\{RG\\}_02_deduped.tsv`
    - | For each readgroup (RG) report the number of unique (deduplicated) reads mapped against the reference genome, the duplication rate
      | and the number of basepairs covered in the reference genome by the reads
      ::

          Order     Family      Species       Reference  ReadsDeduped  DuplicationRate  CoveredBP
          Primates  Hominidae   Homo_sapiens  fixed      98            2.31             4216
  * - :file:`$\\{RG\\}_03_bedfiltered.tsv`
    - | For each readgroup (RG) show the number of reads remaining in the bam-file after bedfiltering as well as the number of covered basepairs
      | in the reference genome
      ::

          Order     Family      Species       Reference  ReadsBedfiltered PostBedCoveredBP
          Primates  Hominidae   Homo_sapiens  fixed      97               4177
  * - :file:`$\\{RG\\}_04_deamination.tsv`
    - For each readgroup (RG) show the deamination stats for the mapped bam-file after bedfiltering
      ::

          Ancientness:  ++  = more than 9.5% of the reads that show a terminal C in both the 5' and 3' position in the reference genome, carry a T
                        +   = more than 9.5% of the reads that show a terminal C in either the 5' or 3' position in the reference genome, carry a T
                        -   = no signs for DNA deamination patterns

          ReadsDeam(1term): The number of reads (after deduplication and bedfiltering) that show a deamination in the terminal base positions
          ReadsDeam(3term): The number of reads (after deduplication and bedfiltering) that show a deamination in the three terminal base positions
          Deam5(95ci):      For the terminal 5' end, the percentage of C to T substitutions (and the 95% confidence interval)
          Deam3(95ci):      For the terminal 3' end, the percentage of C to T substitutions (and the 95% confidence interval)
          Deam5Cond(95ci):  Taken only 3' deaminated sequences, report the percentage of C to T substitutions (and the 95% confidence interval) at the 5' terminal base
          Deam3Cond(95ic):  Taken only 5' deaminated sequences, report the percentage of C to T substitutions (and the 95% confidence interval) at the 3' terminal base

**final report:**

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - File
     - Description

   * - :file:`final_report.tsv`
     - A summary of all the files in the :file:`stats` dir plus additional information gathered from processes during the pipeline run::

        RG                  The analyzed readgroup
        ReadsRaw            The total number of reads in the split' file (by RG)
        ReadsFiltered       The total number of reads in the split' file after the bamfilter process - removing paired reads. (by RG)
        ReadsLengthfiltered The total number of reads in the split' file after the filterLength process - removing reads <35bp length (by RG)
        SpeciesKmers        The total number of kmers assigned to the species node by krakenuniq
        KmerCoverage        The proportion of species kmers assigned by the total number of kmers present in the database for that species
        KmerDupRate         The average duplication rate of kmers used for the assignment of the given species
        ExtractLVL          The taxon level used for extraction of reads (f or o)
        ReadsExtracted      The number of reads extracted/assigned to the taxon by krakenuniq
        Order               The name of the Order assigned
        Family              The name of the Family assigned
        Species             The name of the reference species used for mapping assigned taxon reads against
        Reference           'fixed' if reference genome is set for that family, 'best' if inferred from the kraken assignments
        ReadsMapped         The number of reads mapped against the Species genome
        ProportionMapped    The proportion of mapped/extracted reads
        ReadsDeduped        The number of deduplicated(unique) reads in the alignment file
        DuplicationRate     The duplication rate of mapped reads --> mapped/deduped reads
        CoveredBP           The number of basepairs in the reference genome covered by the aligned reads
        ReadsBedfiltered    The number of reads not overlapping low-complexity regions
        PostBedCoveredBP    The number of basepairs in the reference genome covered by the aligned reads - after bedfiltering
        FamPercentage       Taken all bedfiltered reads of the RG, report the percentage of bedfiltered reads for the given family
        Ancientness         One of ++,+ or - --> See above for an explanation of the symbols
        ReadsDeam(1term)    The number of reads showing a C to T substitution on either of the 5' or 3' ends in respect to the reference
        ReadsDeam(3term)    The number of reads showing a C to T substitution in the terminal 3 basepairs in respect to the reference
        Deam5(95ci)         For the given family, the percentage of C to T substitutions (and the 95% confidence interval) on the terminal 5' end
        Deam3(95ci)         For the given family, the percentage of C to T substitutions (and the 95% confidence interval) on the terminal 3' end
        Deam5Cond(95ci)     Taken only 3' deaminated sequences, report the percentage of C to T substitutions (and the 95% confidence interval) at the 5' terminal base
        Deam3Cond(95ci)     Taken only 5' deaminated sequences, report the percentage of C to T substitutions (and the 95% confidence interval) at the 3' terminal base

| The :file:`cc_estimates.tsv` files contains information about index-hopping and cross contamintaion
| The :file:`nextflow` directory contains nextflow specific information about the run
| the :file:`work` directory can be deleted after the run - it contains nextflow specific intermediate files
