.. _usage-page:

Usage
=====

Flags
-----
Required
""""""""

| **Database-Flags:**
| (see the :ref:`set-up section <setup>` to create these datastructures)

--db DIR    The preindexed ``Kraken`` or ``KrakenUniq`` database used for the classification of reads
--genome DIR    The directory containing the reference genomes in the format :file:`DIR/$family/$species.fasta`
--bedfiles DIR  The directory containing the dustmasked bedfiles in the format :file:`DIR/$species_masked.bed`

| **Input-Flags:**
| (See :ref:`input`)

--split DIR     The directory containing the input files in the format :file:`DIR/*.\\{bam|fastq\\}`
--bam BAM-FILE  A still multiplexed bam-file, containing all the reads
--rg TSV-FILE   The expected readgroups and primer combinations to demultiplex the BAM_FILE by. Fileformat: :code:`RG\\tP7\\p5`

.. note::
    
    For the the internal use (MPI EVA), there is a demultiplexer available in the pipeline. It is thus possible to run quicksand with :code:`-bam BAM-FILE` and :code:`--rg TSV-FILE` as alternative to the :code:`--split` parameter. 


Optional
""""""""
**Pipeline parameters**

--analyze   Given this flag, quicksand will check all alignments in :file:`bedfilteredReads.bam` for C to T substitutions compared to the reference genome. The file :file:`$RG_deamination_stats.tsv` is created for each readgroup. See the :ref:`files` section. This step takes some time.
--report    Provided with this flag, quicksand will output a 'final' summary called :file:`final_report.tsv` that summarizes all the stats scattered across the :file:`stats` directory. 
--byrg  Change the order of output-files in the :file:`out` directory from a taxonomy-based structure to a readgroup-based structure. See the :ref:`output` section for a comparison  
--specmap TSV-FILE   With the TSV_FILE containing families and species in the format ``Family\tSpecies_name,Species_name`` in each line, ignore the "best taxon" assigned by the kraken classifier as reference for mapping for given families. Instead map the :file:`extractedReads` of the $family to the species specified in the TSV-FILE. The ``Species_name`` must correspond to a filename (without :code:`.fasta`) in the :file:`genomes/$family/` directory
--capture STRING    Given a STRING ``Family,Family,...``, apply reduced filters to the given families. For provided families, the bedfiltering is ommited


**Process parameters**

--cutoff N  Length cutoff after BWA-mapping (default: 35). Remove all reads below a length of N from the alignemnt.
--quality N     Mapping quality filter after BWA-mapping (default:25). Remove all reads with a mapping-quality of <N from the alignment.
--min_kmers N   Minimum required unique kmers assigned by KrakenUniq (default: 129). To reduce false-positives: Remove families from the assignments if the total number of assigned unique kmers to taxa of that family by KrakenUniq is below N
--min_reads N   Minimum required assigned reads by KrakenUniq (default: 3). To reduce false-positives: Remove families from the assignments if the total number of assigned reads to taxa of that family by KrakenUniq is below N
--keeppaired    If :file:`.bam` files are provided as input, by default, unmerged reads are filtered out. This flag prevents that.
--filterunmapped    If premapped :file:`.bam` files are provided as input, use this flag to filter out unmapped reads from the bam-files
--level N   Set BGZF compression level (default: 6)
--krakenthreads N   Number of threads per Kraken process (default: 4)

Nextflow
""""""""

A selection of built-in Nextflow flags that may be of use (Be aware: only one dash - with these flags)::

    -profile  <profile>  pick a profile for executing the pipeline (see below)
    -resume              Resume processing; do not re-run completed processes
    -N        <email>    send notification to email upon fiishing
    -c        <file>     path to additional nextflow.config file for local parameters
    -w        <path>     specify the "work" directory for the nextflow intermediate files


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
    Quicksand will process all the :file:`.fastq` and :file:`.bam`-files within the :code:`--split` directory and ignore all remaining files.
    Be sure to name/rename your files accordingly (e.g rename :file:`.fq` to :file:`.fastq` files)


Run the pipeline
--------------------

As outlined below, there are several ways to run the pipeline. 
Please see to the :ref:`installation <install-page>` to set up the required underying databases and 
the :ref:`input section <input>` to know about the required input.

.. attention::

    Before starting the pipeline make sure that you are in separate folder. The pipeline writes its output into the current working directory! 


Run from local repository
"""""""""""""""""""""""""

With everything set up quicksand can be executed by pointing nextflow to the :file:`main.nf` file in the cloned repository::

    mkdir runDir && cd runDir
    export NXF_SINGULARITY_CACHEDIR=~/pipeline/singularity/ 

    nextflow run ~/pipeline/quicksand/main.nf \
        --split     PATH/TO/INPUT/DATA \
        --genome    ~/pipeline/data/genomes \
        --db        ~/pipeline/data/kraken/Mito_db_kmer22 \
        --bedfiles  ~/pipeline/data/masked \
        -profile    singularity
        --report 
        --analyze 


Run from remote repository
""""""""""""""""""""""""""
Instead of cloning the repository, it is also possible to run the pipeline directly from github::

    mkdir runDir && cd runDir
    export NXF_SINGULARITY_CACHEDIR=~/pipeline/singularity/ 

    nextflow run mpieva/quicksand \
        --split     PATH/TO/INPUT/DATA \
        --genome    ~/pipeline/data/genomes \
        --db        ~/pipeline/data/kraken/Mito_db_kmer22 \
        --bedfiles  ~/pipeline/data/masked \
        -profile    singularity
        --report 
        --analyze


It is possible to specify a version (or branch) by providing the :code:`-r <branch/tag>` flag. By default the :code:`master` branch is pulled 
and stored locally in the hidden :file:`~/.nextflow` directory. If the remote pipeline is updated, make sure to pull the updates by running::

    nextflow pull mpieva/quicksand

To reduce the number of flags manually handed over to the pipeline, see the :ref:`configuration section <configuration-page>` on how to set up an individual configuration of the pipeline.

.. _output:

Output
------

For each readgroup, the pipeline outputs the raw reads in :file:`.bam`-format at each stage of the pipeline. 
Additionally the number of reads assigned to one family/extracted, mapped, deduplicated, and bedfiltered as well as the estimated DNA-damage
are reported for a quick overview in one big summary-file.

Two ways exist to structure the output. Based on the :code:`--byrg` flag either by readgroup or by family.
without the :code:`--byrg` flag, the output is structured by Family

Structured by Family
""""""""""""""""""""

This overview corresponds to a run with the :code:`--analyze` and the :code:`--report` flags provided.
Several directories and files should appear after the run::

    RunDir
    ├── out
    │    └── {family}
    │         ├── {readgroup}_extractedReads-{family}.bam
    │         ├── aligned
    │         │    ├── {readgroup}.{species}.bam
    │         │    └── {readgroup}.{species}_deduped.bam
    │         └── bed
    │              └── {readgroup}.{species}_deduped_bedfiltered.bam
    ├── kraken
    │    ├── {readgroup}.report
    │    └── {readgroup}.translate
    ├── stats
    │    ├── splitcounts.tsv
    │    ├── {readgroup}_extracted.tsv
    │    ├── {readgroup}_mapped.tsv
    │    ├── {readgroup}_mapped_coverage.tsv
    │    ├── {readgroup}_unique_mapped.tsv
    │    ├── {readgroup}_bedfiltered.tsv
    │    └── {readgroup}_deamination_stats
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
    │         ├── {readgroup}_extractedReads-{family}.bam
    │         ├── aligned
    │         │    ├── {family}.{species}.bam
    │         │    └── {family}.{species}_deduped.bam
    │         └── bed
    │              └── {family}.{species}_deduped_bedfiltered.bam
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
   * - :file:`extractedReads.bam`
     - For the given readgroup and family, all reads assigned by KrakenUniq to that family
   * - :file:`alignedReads.bam`
     - For a given family and species, all extractedReads (see above) mapped/aligned to the genome of the assigned species
   * - :file:`aligned_dedupedReads.bam`
     - The aligned reads, but depleted of PCR duplicates
   * - :file:`bedfilteredReads`
     - The aligned_deduped reads but additionally depleted of reads overlapping low-complexity regions


**kraken/**

.. list-table::
   :widths: 20 80
   :header-rows: 1
   
   * - File
     - Description
   * - :file:`kraken.report`
     - The raw KrakenUniq report
   * - :file:`kraken.translate`
     - The human readable kraken report in mpa-format

**stats/**

.. list-table::
   :widths: 20 80
   :header-rows: 1
   
   * - File
     - Description
   * - :file:`splitcounts.tsv`
     - Contains the number of reads per readgroup
   * - :file:`extracted.tsv`
     - For each readgroup show the number of reads per extracted family
   * - :file:`mapped.tsv`
     - Show the number of reads mapped against the reference genome
   * - :file:`mapped_coverage.tsv`
     - Contains the number of covered basepairs for each mapping
   * - :file:`unique_mapped.tsv`
     - Contains the number of deduplicated reads mapped against the reference genome
   * - :file:`bedfiltered.tsv`
     - The number of reads remaining in the bam-file after bedfiltering
   * - :file:`deamination_stats`
     - Show an estimation of the 'ancientness' of families based on deamination frequencies of recovered reads
     
| The :file:`final_report.tsv` is a summary of all the files in the :file:`stats` dir.
| The :file:`report` directory contains information about the run
| the :file:`work` directory can be deleted after the run - it contains nextflow specific intermediate files 
