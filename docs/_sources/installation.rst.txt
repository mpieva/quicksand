.. role:: bold
.. role:: heading1
.. _install-page:

Installation
============

Use this section to set up the pipeline to work with real data. To test if the pipeline works on your computer, see the :ref:`quickstart-page` section

.. _requirements:

Requirements
------------

Make sure you have the following two components installed.

:Nextflow: See more details here: `Nextflow <https://www.nextflow.io/docs/latest/getstarted.html>`_
:Singularity or Docker: See more details here: `Singularity <https://sylabs.io/guides/3.0/user-guide/installation.html>`_ or `Docker <https://docs.docker.com/get-docker/>`_

Nextflow is the pipeline framework, while Singularity/Docker are software tools for the containerization of code. 

.. tip::
    
    check the successful installation of the software by running::
        
        nextflow -v
            >>> nextflow version 20.04.1.5335
        singularity --version
            >>> singularity version 3.7.2-dirty


:heading1:`Working with test data`

Test runability
---------------

To test the runability of the pipeline, see the :ref:`quickstart-page` section

.. _container:

Container Cache
---------------

Quicksand uses container software (Docker,Singularity) to ensure the pipeline runs within a stable environment. 
For single-tool processes, container images are either pulled from the `Galaxy image repository <https://depot.galaxyproject.org/singularity>`_ or the
`Quai.io biocontainer repository <https://quay.io/organization/biocontainers>`_. For multi-tool processes and
custom functions, self-built images are hosted on `Dockerhub <https://hub.docker.com/r/merszym>`_. The 
underlying Dockerfiles for those images can be found within the :file:`assets/docker` directory of the `repository <https://www.github.com/mpieva/quicksand>`_.

To reuse downloaded images for multiple pipeline runs, specify a cache-directory. type::

    mkdir singularity
    export NXF_SINGULARITY_CACHEDIR=$PWD/singularity

.. note::
   Make sure that Singularity has the permissions to overlay and access your file-system within the 
   container. Otherwise the pipeline wont be able to read from or write to it. Upon a 
   :bold:`"file doesn't exist" error`, create an additional :file:`nextflow.config` config-file

   Add the following content to the file::
    
        singularity {
          runOptions = "--bind /directory/in/use"
        }
   
   Add this file to your nextflow run with the :code:`-c` flag::
    
        nextflow run ... -c singularity/nextflow.config

.. attention::
    the :code:`-profile` and the :code:`-c` flag have only one dash!

Now -again- the pipeline can be tested by running::

    nextflow run -profile test,singularity -c singularity/nextflow.config

| The meaning of the flags and the different ways of customizing the pipeline is described in the :ref:`usage-page` section. 
| In case of choosing Docker over Singularity, use the :code:`-profile test,docker` command.

The :file:`singularity` directory should now contain all the images used in the pipeline::

    singularity
    ├── depot.galaxyproject.org-singularity-samtools-1.15.1--h1170115_0.img
    ├── merszym-biohazard_bamrmdup-v0.2.img
    └── merszym-quicksand-1.2.img

:heading1:`Working with real data`

.. _setup:

Create datastructure
--------------------

The required underlying datastructure of the pipeline is in detail described in the :ref:`quicksand_build-page` section

In short: You need a :bold:`precompiled kraken database`, the respective :bold:`reference genomes` and :bold:`bedfiles` indicating low-complexity regions. 
Use the supplementary pipeline :code:`quicksand-build` (once) to download the taxonomy from NCBI/taxonomy, all mitochondrial
genomes from NCBI/RefSeq and create the required databases and files for you.

For this session create the datastructure for the :bold:`Primate mtDNA` from RefSeq::
	
	nextflow run mpieva/quicksand-build --outdir refseq --include Primates

.. attention::
    | Building the database requires ~40G of RAM
    | Be patient, downloading the taxonomy plus the creation of the database might take :bold:`~1h`.


This command creates a directory :file:`refseq` that contains the files required to run quicksand::

    refseq
    ├── kraken
    │    └── Mito_db_kmer22
    ├── genomes
    │    ├── {family}
    │    │    └── {species}.fasta
    │    └── taxid_map.tsv
    └── masked
         └── {species}_masked.bed


With the datastructure created, the pipeline is ready to be used with the following flags::

    --db         refseq/kraken/Mito_db_kmer22/
    --genomes    refseq/genomes/
    --bedfiles   refseq/masked/   



Run the pipeline
----------------

As :bold:`input` for the pipeline, download the Hominin "Hohlenstein-Stadel" mtDNA [1]_ into a directory :bold:`split`
::
	
	wget -P split http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/BAM/mtDNA/HST.raw_data.ALL.bam

And run the quicksand pipeline::

    nextflow run mpieva/quicksand \
        --db        refseq/kraken/Mito_db_kmer22/ \
        --genomes   refseq/genomes/ \
        --bedfiles  refseq/masked/ \
        --split     split \
        -profile    singularity

| Please see the :ref:`usage-page` section for an explaination of the flags and the input!
| Please see the :ref:`output` section for an explaination of the output files!

A summary of all the stats can be found in the :file:`final_report.tsv` file


Filter the Results
------------------

As can be seen in the :code:`final_report.tsv`, not all sequences were assigned to Homindae, but to a couple of other Primate families too.
The assignment of false positive taxa is a well-known problem of kmer-based assignment methods and additional filters need to be applied.

Based on simulated data, our recommended cutoffs are:

- :bold:`FamPercentage` cutoff of 1% and/or
- :bold:`ProportionMapped` cutoff of 0.5-0.7.

The kmer-information is also indicative. If the :bold:`FamilyKmers` and :bold:`KmerCoverage` values are low and
the :bold:`KmerDupRate` value is high, the assigment of the family is only based on a small number of kmers within the reads 


.. [1] http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/README
