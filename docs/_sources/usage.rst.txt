.. role:: bold
.. role:: mono

.. _usage-page:

Usage
=====

This section explains the input and output of the pipeline, as well as the commands and flags to run it

Required Flags
--------------

| **Datastructure flags**

.. list-table::
  :widths: 10 10 60
  :header-rows: 1
   
  * - Flag
    - Input type
    - Description
  * - :mono:`--db`
    - PATH
    - The :bold:`directory` containing the preindexed :file:`Kraken` or :file:`KrakenUniq` database (see: :ref:`quicksand_build-page`):: 

        Input:

        refseq
          ├── kraken
          │    └── Mito_db_kmer22
          │           ├── taxonomy
          │           ├── ...
          │           └── database.kdb  


        Example:

        --db Path/to/refseq/kraken/Mito_db_kmer22 

  * - :mono:`--genomes`
    - PATH
    - | The :bold:`directory` containing the indexed FASTA-FILES of the reference genomes that were used to build
      | the kraken database. Format :file:`PATH/$\\{family\\}/$\\{species\\}.fasta` (see: :ref:`quicksand_build-page`)
      ::

          Input:

          refseq
            ├── genomes
            │    └── Hominidae
            │           ├── Homo_sapiens.fasta
            │           ├── Homo_sapiens.fasta.fai
            │           └── ...


          Example:

          --genomes Path/to/refseq/genomes

  * - :mono:`--bedfiles`
    - PATH
    - | The :bold:`directory` containing the dustmasked BED-FILES of the reference genomes FASTA-FILES
      | Format :file:`DIR/$\\{species\\}_masked.bed` (see: :ref:`quicksand_build-page`)
      ::

          Input:

          refseq
            ├── masked
            │      ├── Homo_sapiens_masked.bed
            │      └── ...


          Example:

          --bedfiles Path/to/refseq/masked


| **Input flags**

.. note::
    
  The pipeline contains an optional initial demultiplexing process. However, that demultiplexer works only on bam-files 
  provided by the `MPI EVA Core Unit <https://www.eva.mpg.de/de/genetics/index/>`_ as the result of library preparation and sequencing.
  For the processing of data coming from the MPI EVA run quicksand with the :code:`--bam PATH` and :code:`--rg PATH` flags as alternative 
  to the :code:`--split PATH` parameter. 


.. list-table::
  :widths: 10 10 60
  :header-rows: 1
   
  * - Flag
    - Input type
    - Description
  * - :mono:`--split`
    - PATH
    - | Standard input
      | A **directory** containing the demultiplexed, adapter trimmed (and overlap-merged) input files 
      ::

          Input:

          split
            ├── RG1.bam
            ├── RG2.fastq
            ├── RG3.fastq.gz
            ├── RG4.fq.gz
            └── ...

          Example:

          --split Path/to/split/

  * - :mono:`--bam`
    - PATH
    - | Use together with the :code:`--rg` flag
      | The multiplexed BAM-FILE, as provided by the `MPI EVA Core Unit <https://www.eva.mpg.de/de/genetics/index/>`_ containing 
      | adapter-trimmed and overlap-merged sequencing reads
      ::

          Example:

          --bam Path/to/input.bam

  * - :mono:`--rg`
    - PATH
    - | Use together with the :code:`--bam` flag
      | A TSV-FILE, containing library information for the demultiplexing step.
      | Provide the readgroups and respective primer combinations contained in the BAM FILE
      ::

          Input (index.tsv):
          
          #Index library ID	primer_P7	primer_P5
          RG1	1113	1137
          RG2	1114	1138

          Example:

          --rg Path/to/index.tsv


Optional Flags
--------------

.. list-table::
  :widths: 10 10 60
  :header-rows: 1
  
  * - Flag
    - Input type
    - Description

  * - :mono:`--specmap`
    - PATH
    - | Provide a config TSV file (No header)
      | Map Kraken-assigned family reads (:file:`extractedReads`) to the specified species instead of the 'best taxa' found by the pipeline.
      | Species name must correspond to a filename in :file:`genomes`.
      | Separate multiple species by comma       
      ::
        
          Input (specmap.tsv):
          
          Hominidae Homo_sapiens

          Example:

          --specmap Path/to/specmap.tsv

  * - :mono:`--byrg`
    - 
    - | If provided, change the ordering of the output-files in the :file:`out` directory from a taxonomy-based structure 
      | to a readgroup-based structure. 
      | See the :ref:`output` section for a comparison      
      ::

          Example:

          --byrg

  * - :mono:`--taxlvl`
    - [o,f]
    - | Default: f
      | Bin assigned reads by the specified taxon level (family or order level): e.g. Homindae or Primates.
      | These binned reads (:file:`extractedReads`) are mapped against the genome(s) of each families 'best taxa'.
      | or :code:`--specmap` taxa. e.g. Map all reads assigned to Primates to the Homo_sapiens genome 
      | **Note:** For the order-level bins, :file:`extractedReads` are mapped several times to different (family) genomes.    
      ::

          Example:

          --taxlvl o

  * - :mono:`--skip_analyze`
    - 
    - | Skip the analyzeDeamination process during the pipeline run
      ::

          Example:

          --skip_analyze

  * - :mono:`--skip_report`
    - 
    - | Skip the creation of the :file:`final_report.tsv` during the pipeline run
      ::

          Example:

          --skip_report

  * - :mono:`--skip_bed`
    - STRING
    - | Provide a string of comma-separated families
      | Skip the runIntersectBed process during the pipeline run for the specified families
      ::

          Example:

          --skip_bed Hominidae,Hyaenidae

**Process parameters**

.. list-table::
  :widths: 10 10 60
  :header-rows: 1
  
  * - Flag
    - Input type
    - Description

  * - :mono:`--bamfilterflag`
    - N
    - | For process filterBam
      | Filter the BAM file based on the provided SAMTOOLS FLAG (default: 1 = filter paired reads).
      | see `HERE <https://broadinstitute.github.io/picard/explain-flags.html>`_ to find a desired filterflag
      ::

          Example:

          --bamfilterflag 5
  
  * - :mono:`--bamfilter_length_cutoff`
    - N
    - | For process filterLength
      | Filter out reads below the given length cutoff (default: 35).
      ::

          Example:

          --bamfilter_length_cutoff 35

  * - :mono:`--bamfilter_quality_cutoff`
    - N
    - | For process mapBwa
      | Filter out reads with a mapping quality below the given quality cutoff (default: 25).
      ::

          Example:

          --bamfilter_quality_cutoff 25

  * - :mono:`--krakenuniq_min_kmers`
    - N
    - | For process runKrakenUniq
      | A filter to reduce the number of false-positive krakenuniq assignments
      | Filter out assigned families from the krakenuniq classification with a kmer-count < N (default: 129).
      ::

          Example:

          --krakenuniq_min_kmers 129

  * - :mono:`--krakenuniq_min_reads`
    - N
    - | For process runKrakenUniq
      | A filter to reduce the number of false-positive krakenuniq assignments
      | Filter out families from the krakenuniq classification with < N reads assigned (default: 3).
      ::

          Example:

          --krakenuniq_min_reads 3

  * - :mono:`--compression_level`
    - [0-9]
    - | For the BAM output processes
      | Set BGZF compression level (default: 0)
      ::

          Example:

          --compression_level 9


Profiles
--------

| Quicksand includes several profiles that can be used with the :code:`-profile` flag (Be aware: only one dash -)
| delimit multiple profiles by comma

.. list-table::
  :widths: 20 60
  :header-rows: 1
  
  * - Profile
    - Description

  * - singularity
    - Use Singularity as container software

  * - docker
    - Use Docker as container software

  * - conda
    - Use Conda environments for prcesses (WIP)

  * - test
    - Use testdata to run the pipeline :code:`nextflow run mpieva/quicksand -profile test,singularity`

  * - debug
    - Keep intermediate files in the :file:`work` directory

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


Run the pipeline
----------------

There are several ways to run the pipeline. 
Please see to the :ref:`installation <install-page>` to set up the required underying databases and 
the :ref:`input section <input>` to know about the required input.

.. attention::

    The pipeline writes its output into the current working directory! 


Run from local repository
"""""""""""""""""""""""""

To run the pipeline locally, pull the code from github::

  git clone git@github.com:mpieva/quicksand

and executed the pipeline by pointing nextflow to the :file:`main.nf` file in the cloned repository::

    nextflow run quicksand/main.nf \
        --split      /path/to/split/ \
        --genomes    /path/to/refseq/genomes \
        --db         /path/to/refseq/kraken/Mito_db_kmer22 \
        --bedfiles   /path/to/refseq/masked \
        -profile     singularity


Run from remote repository
""""""""""""""""""""""""""

Instead of cloning the repository, run the pipeline directly from github::

    nextflow run quicksand/main.nf [-r v1.5] \
        --split      /path/to/split/ \
        --genomes    /path/to/refseq/genomes \
        --db         /path/to/refseq/kraken/Mito_db_kmer22 \
        --bedfiles   /path/to/refseq/masked \
        -profile     singularity


Specify a version (or branch) by providing the :code:`-r <branch/tag>` flag. By default the :code:`master` branch is pulled 
and stored locally in the hidden :file:`~/.nextflow/assets/mpieva/quicksand` directory. If the remote pipeline is updated, make sure to pull the updates by 
running::

    nextflow pull mpieva/quicksand 

To reduce the number of flags manually handed over to the pipeline, see the :ref:`configuration section <configuration-page>` on how to set up an individual configuration of the pipeline.

.. _output:

Output
------

Based on the :code:`--byrg` flag, structure the output by either readgroup or family. By default, the output is structured by Taxon


Structured by Taxon
""""""""""""""""""""

Several directories and files should appear after the run. The 'taxon' corresponds to either the family or the order level name::

    RunDir
    ├── out
    │    └── {taxon}
    │         ├── {readgroup}_extractedReads-{taxon}.bam
    │         ├── aligned
    │         │    ├── {family}.{taxon}.{species}.bam
    │         │    └── {family}.{taxon}.{species}_deduped.bam
    │         ├── bed
    │         │    └── {family}.{taxon}.{species}_deduped_bedfiltered.bam
    │         └── deaminated
    │              └── {family}.{taxon}.{species}_deduped_deaminated53.bam
    ├── kraken
    │    ├── {readgroup}.report
    │    └── {readgroup}.translate
    ├── stats
    │    ├── splitcounts.tsv
    │    ├── {readgroup}_extracted.tsv
    │    ├── {readgroup}_mapped.tsv
    │    ├── {readgroup}_mapped_coverage.tsv
    │    ├── {readgroup}_mapped_deduped.tsv
    │    ├── {readgroup}_mapped_deduped_bedfiltered.tsv
    │    └── {readgroup}_mapped_deduped_deamination53.tsv
    ├── reports
    │    ├── report.html
    │    ├── timeline.html
    │    └── trace.tsv
    ├── work
    │    └── ...
    └── final_report.tsv 

Structured by Readgroup
"""""""""""""""""""""""

This overview of the :file:`out/` dir corresponds to a :code:`--byrg` run::

    RunDir
    ├── out
    │    └── {readgroup}
    │         ├── {readgroup}_extractedReads-{taxon}.bam
    │         ├── aligned
    │         │    ├── {readgroup}.{family}.{species}.bam
    │         │    └── {readgroup}.{family}.{species}_deduped.bam
    │         ├── bed
    │         │    └── {readgroup}.{family}.{species}_deduped_bedfiltered.bam
    │         └── deaminated
    │              └── {readgroup}.{family}.{species}_deduped_deaminated53.bam
    ...

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
   * - :file:`deaminated/$\\{RG\\}.$\\{family\\}.$\\{species\\}_deduped_deaminated53.bam`
     - | The same aligned, deduped and bedfiltered bamfile, filtered for reads that show a C to T 
       | substitution at one of the terminal basepair positions in respect to the reference genome
 


**kraken/**

.. list-table::
   :widths: 20 80
   :header-rows: 1
   
   * - File
     - Description
   * - :file:`$\\{RG\\}.report`
     - The standard krakenuniq report
   * - :file:`$\\{RG\\}.translate`
     - The read-wise human readable kraken report in mpa-format

**stats/**

.. list-table::
   :widths: 20 80
   :header-rows: 1
   
   * - File
     - Description
   * - :file:`stats/splitcounts.tsv`
     - For each readgroup (RG), show the number of reads before and after the filterBam process::

        readgroup   split count   filtered count
        test1       235           235
        test2       235           235
        test3       235           235

   * - :file:`$\\{RG\\}_extracted.tsv`
     - For each readgroup (RG) show the number of reads extracted for each assigned taxon::

        Taxon       ReadsExtracted
        Hominidae   235

   * - :file:`$\\{RG\\}_mapped.tsv`
     - For each readgroup (RG) show the number of reads mapped against the reference genome::

        Order     Family      Species       ReadsMapped
        Primates  Hominidae   Homo_sapiens  235
      
   * - :file:`$\\{RG\\}_mapped_coverage.tsv`
     - For each readgroup (RG) show the number of covered basepairs in the reference for each mapping::

        Order     Family      Species       CoveredBP
        Primates  Hominidae   Homo_sapiens  4616        

   * - :file:`$\\{RG\\}_mapped_deduped.tsv`
     - For each readgroup (RG) report the number of unique (deduplicated) reads mapped against the reference genome::

        Order     Family      Species       ReadsDeduped
        Primates  Hominidae   Homo_sapiens  98  

   * - :file:`$\\{RG\\}_mapped_deduped_bedfiltered.tsv`
     - For each readgroup (RG) show the number of reads remaining in the bam-file after bedfiltering::

        Order     Family      Species       ReadsBedfiltered
        Primates  Hominidae   Homo_sapiens  97

   * - :file:`$\\{RG\\}_mapped_deduped_deaminated53.tsv`
     - For each readgroup (RG) show the deamination stats for the mapped bam-file after bedfiltering::

        Ancient:  ++  = more than 9.5% of the reads that show a terminal C in both the 5' and 3' position in the reference genome, carry a T
                  +   = more than 9.5% of the reads that show a terminal C in either the 5' or 3' position in the reference genome, carry a T
                  -   = no signs for DNA deamination patterns
        
        Deam5(95CI): For the terminal 5' end, the percentage of C to T substitutions (and the 95% confidence interval) 
        Deam3(95CI): For the terminal 3' end, the percentage of C to T substitutions (and the 95% confidence interval) 
        Deam5Cond(95CI): Taken only 3' deaminated sequences, report the percentage of C to T substitutions (and the 95% confidence interval) at the 5' terminal base
        Deam3Cond(95CI): Taken only 5' deaminated sequences, report the percentage of C to T substitutions (and the 95% confidence interval) at the 3' terminal base

        Order       Family      Species         Ancient   ReadsDeaminated   Deam5(95CI)           Deam3(95CI)         Deam5Cond(95CI)     Deam3Cond(95CI)
        Primates    Hominidae   Homo_sapiens    ++        13                28.57(13.22,48.67)    27.27(10.73,50.22)  50.0(1.26,98.74)    50.0(1.26,98.74)

**final report:**

.. list-table::
   :widths: 20 80
   :header-rows: 1
   
   * - File
     - Description

   * - :file:`final_report.tsv`
     - A summary of all the files in the :file:`stats` dir plus additional information gathered from processes during the pipeline run::

        RG                  The analyzed readgroup 
        FamilyKmers         The total number of kmers assigned to the family node by krakenuniq
        KmerCoverage        The proportion of family kmers assigned by the number of kmers present in the database for that family
        KmerDupRate         The average duplication rate of kmers used for the assignment of the given family  
        Order               The name of the Order assigned
        Family              The name of the Family assigned
        Species             The name of the reference species used for mapping assigned taxon reads against 
        ExtractLVL          The taxon level used for extraction of reads (f or o)
        ReadsExtracted      The number of reads extracted/assigned to the taxon by krakenuniq
        ReadsMapped         The number of reads mapped against the Species genome  
        ProportionMapped    The proportion of mapped/extracted reads
        ReadsDeduped        The number of deduplicated(unique) reads in the alignment file 
        DuplicationRate     The duplication rate of mapped reads --> mapped/deduped reads
        CoveredBP           The number of basepairs in the reference genome covered by the aligned reads
        ReadsBedfiltered    The number of reads not overlapping low-complexity regions
        FamPercentage       Taken all bedfiltered reads of the RG, report the percentage of bedfiltered reads for the given family
        Ancientness         One of ++,+ or - --> See above for an explanation of the symbols
        ReadsDeaminated53   The number of reads showing a C to T substitution on either of the 5' or 3' ends in respect to the reference
        Deam5               For the given family, the percentage of C to T substitutions (and the 95% confidence interval) on the terminal 5' end
        Deam3               For the given family, the percentage of C to T substitutions (and the 95% confidence interval) on the terminal 3' end
        Deam5Cond           Taken only 3' deaminated sequences, report the percentage of C to T substitutions (and the 95% confidence interval) at the 5' terminal base
        Deam3Cond           Taken only 5' deaminated sequences, report the percentage of C to T substitutions (and the 95% confidence interval) at the 3' terminal base


| The :file:`report` directory contains nextflow specific information about the run
| the :file:`work` directory can be deleted after the run - it contains nextflow specific intermediate files 
