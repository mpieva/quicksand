Installation
============

This section walks you through the stages of downloading, installing and
testing the pipeline. Please follow the instructions step by step.

01. Install Requirements
------------------------

Before coming to that, make sure you have the two following programs installed.

:Nextflow: See more details here: `Nextflow <https://www.nextflow.io/docs/latest/getstarted.html>`_
:Singularity: See more details here: `Singularity <https://sylabs.io/guides/3.0/user-guide/installation.html>`_

    or Docker: See `here <https://docs.docker.com/get-docker/>`_

Nextflow is the language the pipeline is written in, while Singularity/Docker are software tools for
containerization of code - used to run software in a consistent environment. 

.. tip::
    
    check the installation of the software above by running::
        
        nextflow -v
            >>> nextflow version 20.04.1.5335
        singularity --version
            >>> singularity version 3.7.2-dirty


02. Download Repository
-----------------------

The code for the sediment_nf pipeline is hosted on github. The repository can be downloaded
directly from there. The code can be found `here <https://github.com/MerlinSzymanski/sediment_nf/>`_.

Please open your terminal and type::
    
    cd
    mkdir pipeline && cd pipeline
    git clone https://www.github.com/MerlinSzymanski/sediment_nf
    
This code creates a new directory :file:`pipeline` that will contain
all the data produced over the course of this setup process. 
Here - the cloned github repository.


03. Test Runabilty
------------------

With the requirements met, the pipeline can be tested.
To do that, unpack the test database provided in the repository::

    mkdir testrun && cd testrun
    tar -xvzf ../sediment_nf/assets/test/kraken/database.tar.gz

This will create a directory :file:`TestDB` in the current folder. 

To ensure a stable environment, the pipeline runs within a container that gets
downloaded from the internet, from `Dockerhub <https://hub.docker.com/r/merszym/sediment_nf>`_. 
The Dockerfile used to build that image can be found within the :file:`sediment_nf/docker` directory of the repository.
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
    nextflow run ../sediment_nf/main.nf \
        --split ../sediment_nf/assets/test/split/ \
        --genome ../sediment_nf/assets/test/genomes/ \
        --bedfiles ../sediment_nf/assets/test/masked/ \
        --db TestDB/ \
        --specmap ../sediment_nf/assets/test/genomes/specmap.tsv \
        --analyze \
        --report \
        -profile singularity \
        -c ../singularity/nextflow.config

The meaning of the flags and the different ways of customizing the pipeline is described in the customization section. 
In case of choosing Docker over Singularity, exchange :code:`-profile singularity` with :code:`-profile docker`.  

.. attention::
    the :code:`-profile` flag has only one dash!

If the run was successful, several new files and directories will appear in your current working directory. To
see an explanation of the files, see the output section.


04. Setup Datastructure
-----------------------

To run the pipeline with a real databases a certain datastructure is required.

- A preindexed Kraken-database
- All Mammalian mitochondrial reference genomes from RefSeq in a fasta-format
- Bedfiles for these genomes
- A textfile that points to all species of a clade specified by the NCBI taxID

Instead of creating this structure manually, a different pipeline is used
for that

.. seealso::
    Refer to the README of `that pipeline <https://github.com/MerlinSzymanski/datastructure_nf/>`_ for custom
    stettings of the data structure (e.g. kmer-sizes) and a more detailed explanation of the output.

Type into the terminal::

    cd ..
    git clone https://github.com/MerlinSzymanski/datastructure_nf/
    nextflow run datastructure_nf/main.nf -profile singularity --outdir data 

.. attention::

    The creation of the kraken-datasets requires a lot of RAM. 
    If the pipeline fails, make sure the computer fits the requirements!

Upon finishing creates a folder "data" that contains all the files required to run the pipeline against::

    data/
        kraken/
            Mito_db_kmer22/
        genomes/
            {family}/{species}.fasta
            taxid_map.tsv
        masked/
            {species}.masked.bed

This datastructure can be handed over to the pipeline with the following flags::

    --db         /path/to/data/kraken/Mito_db_kmer22/
    --genome     /path/to/data/genomes/
    --bedfiles   /path/to/data/masked/    

With that being done, please see the "Usage" section





