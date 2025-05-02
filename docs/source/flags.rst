.. _flags-page:

Flags
=====
Required flags
--------------

| **Reference Database**

These flags are related to the `quicksand-build <https://github.com/mpieva/quicksand-build>`_ output.

.. list-table::
  :widths: 20 10 60
  :header-rows: 1

  * - Flag
    - Type
    - Description

  * - --db
    - PATH
    - A **directory** containing the preindexed :file:`Kraken` or :file:`KrakenUniq` database::

        Input:

        refseq
          ├── kraken
          │    └── Mito_db_kmer22
          │           ├── taxonomy
          │           ├── ...
          │           └── database.kdb


        Example:

        --db Path/to/refseq/kraken/Mito_db_kmer22

  * - --genomes
    - PATH
    - | A **directory** containing the indexed FASTA-FILES of the reference genomes used to build
      | the KrakenUniq database. Format :file:`genomes/$\\{family\\}/$\\{species\\}.fasta`
      ::

          Input:

          refseq
            ├── genomes
            │    └── Hominidae
            │           ├── Homo_sapiens.fasta
            │           └── ...


          Example:

          --genomes Path/to/refseq/genomes

  * - --bedfiles
    - PATH
    - | A **directory** containing dustmasked BED-FILES for all the reference genomes
      | Format :file:`masked/$\\{species\\}_masked.bed`
      ::

          Input:

          refseq
            ├── masked
            │      ├── Homo_sapiens_masked.bed
            │      └── ...


          Example:

          --bedfiles Path/to/refseq/masked


| **Input flags**

The input for quicksand is a directory with user-supplied files in BAM or FASTQ format. Adapter-trimming, overlap-merging and sequence 
demultiplexing need to be performed by the user prior to running quicksand. However, quicksand also implements a BAM demultiplexing software. 
This software is working for the bam-files provided by the `MPI EVA Core Unit <https://www.eva.mpg.de/de/genetics/index/>`_. 

In these cases, quicksand can be used with the :code:`--bam PATH` and :code:`--rg PATH` flags as alternative to the :code:`--split PATH` parameter.

.. list-table::
  :widths: 20 10 60
  :header-rows: 1

  * - Flag
    - Type
    - Description

  * - --split
    - PATH
    - | Default input method!
      | A **directory** containing demultiplexed, adapter trimmed (and overlap-merged) BAM or FASTQ files
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

  * - --bam
    - PATH
    - | Use in combination with the :code:`--rg` flag
      | A multiplexed BAM-FILE, as provided by the `MPI EVA Core Unit <https://www.eva.mpg.de/de/genetics/index/>`_ containing
      | adapter-trimmed and overlap-merged sequencing reads
      ::

          Example:

          --bam Path/to/input.bam

  * - --rg
    - PATH
    - | Use in combination with the :code:`--bam` flag
      | A TSV-FILE, containing library information for the demultiplexing step.
      | Provide the readgroups and respective primer combinations contained in the BAM FILE
      ::

          Input (index.tsv):

          #Index library ID	primer_P7	primer_P5
          RG1	1113	1137
          RG2	1114	1138

          Example:

          --rg Path/to/index.tsv


Optional flags
--------------

.. list-table::
  :widths: 20 10 60
  :header-rows: 1

  * - Flag
    - Type
    - Description

  * - --fixed
    - PATH
    - | Provide a TSV file
      | Map :file:`extractedReads` (binned sequences) of detected families to the reference genomes listed, instead of the ones determinded by quicksand.
      | The 'Tag' is used as 'Species' name in the final summary reports and filenames.
      ::

          Input (fixed.tsv):

          Family    Species/Tag  Genome
          Hominidae Homo_sapiens  /path/to/seq.fa
          Hominidae Other_cool_hominin  /path/to/other_cool_hominin.fa

          Example:

          --fixed Path/to/fixed.tsv

  * - --fixed_bedfiltering
    - -
    - | Use in combination with the `--fixed` flag
      | Set this flag to run dustmasking and bedfiltering for the `--fixed` references as well. 
      | 
      ::

          Example:

          --fixed Path/to/fixed.tsv --fixed_bedfiltering

  * - --rerun
    - -
    - | Rerun quicksand in an already processed folder
      | Works together with the :code:`--fixed` flag
      | Map already binned reads of families/orders to the reference genomes listed in the
      | :code:`--fixed` TSV file.
      | The analysis records are **added** to the existing final_report file
      ::

          Example:

          --fixed Path/to/fixed.tsv --rerun

  * - --taxlvl
    - [o,f]
    - | Default: f
      | Change the taxonomic level for binning sequences after KrakenUniq classification (family or order level).
      |
      ::

          Example:

          --taxlvl o

  * - --doublestranded
    -
    - | Count C-to-T substitutions at the 5' and G-to-A substitutions at the 3' end sequence alignments,
      | Default: Count C-to-T substitutions on both the 5' and 3' ends.
      ::

          Example:

          --doublestranded


**Process parameters**

.. list-table::
  :widths: 20 10 60
  :header-rows: 1

  * - Flag
    - Type
    - Description

  * - --bamfilterflag
    - N
    - | For initial BAM-file filtering
      | Filter each BAM-file based on the provided SAMTOOLS FLAG (default: 1 = filter paired reads).
      | see `HERE <https://broadinstitute.github.io/picard/explain-flags.html>`_ to find desired filterflags
      ::

          Example:

          --bamfilterflag 5

  * - --bamfilter_length_cutoff
    - N
    - | For initial BAM-file filtering
      | Filter out reads below the given length cutoff (default: 35).
      ::

          Example:

          --bamfilter_length_cutoff 35

  * - --krakenuniq_min_kmers
    - N
    - | For removal of KrakenUniq background-identifications
      | Remove families from the KrakenUniq classification results with an unique kmer-count of less than N (default: 129).
      ::

          Example:

          --krakenuniq_min_kmers 129

  * - --krakenuniq_min_reads
    - N
    - | For removal of KrakenUniq background-identifications
      | Remove families from the KrakenUniq classification results with less than N reads assigned (default: 3).
      ::

          Example:

          --krakenuniq_min_reads 3

  * - --bamfilter_quality_cutoff
    - N
    - | Filter BAM files adter the BWA mapping step
      | Remove mapped sequences with a mapping quality below the given quality cutoff (default: 25).
      ::

          Example:

          --bamfilter_quality_cutoff 25

  * - --reportfilter_percentage
    - N
    - | For the creation of the filtered_report.tsv file
      | Filter family assignments from the final_report.tsv that have a FamPercentage value of less or equal than N percent (default: 0.5).
      ::

          Example:

          --reportfilter_percentage 0.5

  * - --reportfilter_breadth
    - N
    - | For the creation of the filtered_report.tsv file
      | Filter family assignments from the final_report.tsv that have a ProportionExpectedBreadth value of less or equal than N (default: 0.5).
      ::

          Example:

          --reportfilter_breadth 0.8

  * - --compression_level
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

  * - debug
    - dont delete intermediate files in the :file:`work` directory after successful quicksand execution

