.. _install-page:

Installation
============

In the case that the :ref:`quickstart-page` didnt work, this section will walk you through the stages of downloading, installing and
testing the quicksand pipeline in more detail.

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


Download repository
-------------------

| The code for the quicksand pipeline is hosted on `Github <https://github.com/mpieva/quicksand>`_.
| Clone the repository by typing:
::
    
    mkdir pipeline && cd pipeline
    git clone https://www.github.com/mpieva/quicksand
    
This code creates a new directory :file:`pipeline` that will contain all the data produced over the course of this setup process. Here - the cloned github repository.


Test runability
---------------

With the requirements met, the pipeline can be tested.
To do that, unpack the test database provided in the repository::

    mkdir testrun && cd testrun
    tar -xvzf ../quicksand/assets/test/kraken/database.tar.gz

This will create a directory :file:`TestDB` in the current folder. 

To ensure a stable environment, the pipeline runs within a container that gets
downloaded from the internet, from `Dockerhub <https://hub.docker.com/r/merszym/quicksand>`_. 
The Dockerfile used to build that image can be found within the :file:`quicksand/docker` directory of the repository.


To save the downloaded image at an accessible place, type::

    cd ..
    mkdir singularity
    export NXF_SINGULARITY_CACHEDIR=$PWD/singularity

.. note::
   In case some files are later stored or read from a shared drive, make sure 
   that Singularity has the permissions to overlay and access that file-system within the 
   container. The pipeline will otherwise not be able to read from or write to it. 
   To avoid a "file doesn't exist" error, create an additional
   config-file and bind the paths into the container::
   
        nano singularity/nextflow.config
    
   Add the following content to the file::
    
        singularity {
          runOptions = "--bind /path/to/shared/disc"
        }
   
   To save and close nano, press :kbd:`ctrl` + :kbd:`x` and confirm with :kbd:`y`.
   
   Add this file to your nextflow run with the :code:`-c` flag::
    
        nextflow run ... -c singularity/nextflow.config
        

Now the pipeline can be tested by running::

    cd testrun
    nextflow run    ../quicksand/main.nf \
        --split     ../quicksand/assets/test/split/ \
        --genome    ../quicksand/assets/test/genomes/ \
        --bedfiles  ../quicksand/assets/test/masked/ \
        --db        TestDB/ \
        --specmap   ../quicksand/assets/test/genomes/specmap.tsv \
        --analyze   \
        --report    \
        -c          ../singularity/nextflow.config

The meaning of the flags and the different ways of customizing the pipeline is described in the :ref:`usage-page` section. In case of choosing Docker over Singularity, add :code:`-profile docker` to the command.  

.. attention::
    the :code:`-profile` and the :code:`-c` flag has only one dash!

If the run was successful, several new files and directories will appear in your current working directory. To see an explanation of the files, see the :ref:`output` section.

.. _setup:

04. Setup Datastructure
-----------------------

To run the pipeline with a real database a certain datastructure is required.

- A preindexed Kraken-database
- All Mammalian mitochondrial reference genomes from RefSeq in a fasta-format
- Bedfiles for these genomes
- A list that points to all species of a clade specified by the NCBI taxID

Instead of creating this structure manually, a different pipeline is used
for that

.. seealso::
    Refer to the README of `that pipeline <https://github.com/mpieva/quicksand-build>`_ for custom
    settings of the data structure (e.g. kmer-sizes) and a more detailed explanation of the output.

The datastructure-pipeline can be started directly from the repository by tying::

    cd ..
    nextflow run mpieva/quicksand-build --outdir data 

.. attention::

    The creation of the preindexed kraken-databases requires ~50GB of RAM. 
    If the pipeline fails, make sure the computer fits the requirements!

This creates a folder "data" that contains all the database files required to run quicksand::

    data
    ├── kraken
    │    └── Mito_db_kmer22
    ├── genomes
    │    ├── {family}
    │    │    └── {species}.fasta
    │    └── taxid_map.tsv
    └── masked
         └── {species}.masked.bed

This datastructure can be used by quicksand with the following flags::

    --db         /path/to/data/kraken/Mito_db_kmer22/
    --genome     /path/to/data/genomes/
    --bedfiles   /path/to/data/masked/    


05. Run real Data
-----------------

Before running the test, make sure you create a new directory::

    mkdir runDir && cd runDir

For this testrun with real data, download the Hohlenstein-Stadel mtDNA (please see the [README]_) ::
    

    wget -P split http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/BAM/mtDNA/HST.raw_data.ALL.bam

And run the quicksand pipeline::

    nextflow run ~/pipeline/quicksand/main.nf \
        --db        ~/pipeline/data/kraken/Mito_db_kmer22 \
        --genome    ~/pipeline/data/genomes \
        --bedfiles  ~/pipeline/data/masked \
        --split     split \
        --report    \
        --analyze   \
        -c          ~/pipeline/singularity/nextflow.config

| Please see the :ref:`usage-page` section for an explaination of the flags and the input!
| Please see the :ref:`output` section for an explaination of the output files!

