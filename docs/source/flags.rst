.. _flags-page:

Flags
=====
Required Flags
--------------

| **Datastructure flags**

.. list-table::
  :widths: 10 10 60
  :header-rows: 1

  * - Flag
    - Input type
    - Description

  * - --db
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

  * - --genomes
    - PATH
    - | The :bold:`directory` containing the indexed FASTA-FILES of the reference genomes that were used to build
      | the kraken database. Format :file:`PATH/$\\{family\\}/$\\{species\\}.fasta` (see: :ref:`quicksand_build-page`)
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
    - | The **directory** containing the dustmasked BED-FILES of the reference genomes FASTA-FILES
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


quicksand contains an optional demultiplexing preprocessing. However, it is an inhouse demultiplexer working only on bam-files
provided by the `MPI EVA Core Unit <https://www.eva.mpg.de/de/genetics/index/>`_. For the processing of data coming from the
MPI EVA run quicksand with the :code:`--bam PATH` and :code:`--rg PATH` flags as alternative to the :code:`--split PATH` parameter.

.. list-table::
  :widths: 10 10 60
  :header-rows: 1

  * - Flag
    - Input type
    - Description

  * - --split
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

  * - --bam
    - PATH
    - | Use together with the :code:`--rg` flag
      | The multiplexed BAM-FILE, as provided by the `MPI EVA Core Unit <https://www.eva.mpg.de/de/genetics/index/>`_ containing
      | adapter-trimmed and overlap-merged sequencing reads
      ::

          Example:

          --bam Path/to/input.bam

  * - --rg
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

  * - --fixed
    - PATH
    - | Provide a TSV file
      | Map :file:`extractedReads` to the specified genome for given families instead of the one determinded by quicksand.
      | The tag is used as 'Species' name in the reports and the filenames.
      ::

          Input (fixed.tsv):

          Family    Species(tag)  Genome
          Hominidae Homo_sapiens  /path/to/seq.fa

          Example:

          --fixed Path/to/fixed.tsv

  * - --rerun
    - -
    - | Run the pipeline in an already processed folder
      | Works together with the :code:`--fixed` flag
      | Map already extracted reads of families or orders to all the species assigned in the
      | :code:`--fixed` references file.
      | These records are **added** to the final_report file
      ::

          Example:

          --fixed Path/to/fixed.tsv --rerun

  * - --taxlvl
    - [o,f]
    - | Default: f
      | Change taxon level (family or order level) of binned sequences after KrakenUniq.
      | Binned reads are still mapped against the genomes of each `families` reference genome.
      |
      | Example: Map all reads assigned to Primates to the Homo_sapiens genome
      | **Note:** For the order-level bins, Binned reads are mapped several times to different (family) genomes.
      ::

          Example:

          --taxlvl o


**Process parameters**

.. list-table::
  :widths: 10 10 60
  :header-rows: 1

  * - Flag
    - Input type
    - Description

  * - --bamfilterflag
    - N
    - | For initial bam file filtering
      | Filter the file based on the provided SAMTOOLS FLAG (default: 1 = filter paired reads).
      | see `HERE <https://broadinstitute.github.io/picard/explain-flags.html>`_ to find a desired filterflag
      ::

          Example:

          --bamfilterflag 5

  * - --bamfilter_length_cutoff
    - N
    - | For initial bam file filtering
      | Filter out reads below the given length cutoff (default: 35).
      ::

          Example:

          --bamfilter_length_cutoff 35

  * - --krakenuniq_min_kmers
    - N
    - | For metagenomic classification
      | Remove families from the KrakenUniq classification results with a kmer-count of less than N (default: 129).
      ::

          Example:

          --krakenuniq_min_kmers 129

  * - --krakenuniq_min_reads
    - N
    - | For metagenomic classification
      | Remove families from the KrakenUniq classification results with less than N reads assigned (default: 3).
      ::

          Example:

          --krakenuniq_min_reads 3

  * - --bamfilter_quality_cutoff
    - N
    - | For after the mapping step
      | Filter out reads with a mapping quality below the given quality cutoff (default: 25).
      ::

          Example:

          --bamfilter_quality_cutoff 25

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
    - Keep intermediate files in the :file:`work` directory

