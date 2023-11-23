.. _configuration-page:

Configuration
=============

.. _config:

Config-file
-----------

| **create a custom config file**
| Since most of the required flags won't change across runs - like the :code:`--db` or the :code:`--genomes` flag, it is possible to create a custom :file:`nextflow.config` file
| to define these parameters as default values for the pipeline.
| Create a config-file :file:`nextflow.config` as outlined below and hand it over to the pipeline with the :code:`-c` flag:
::

    params {
        db         = "path/to/kraken/Mito_db_kmer22"
        genomes    = "path/to/genomes"
        bedfiles   = "path/to/masked/"
    }

    //add singularity-parameters to the pipeline
    //bind: add shared drives (if necessary)
    //cacheDir: where to download the container

    singularity{
        runOptions = "--bind /mnt/"
        cacheDir   = "path/to/singularity/"
    }

Please see the available configuration options in the `nextflow documentation <https://www.nextflow.io/docs/latest/config.html#scope-singularity>`_

And run the pipeline::

    nextflow run mpieva/quicksand -c nextflow.config --split path/to/split -profile singularity

Environmental variables
-----------------------

The following nextflow specific ENV variables can be set::

    NXF_SINGULARITY_CACHEDIR <path> Where to save the pulled Singularity-images
    NXF_WORK                 <path> Corresponds to the -w flag
    NXF_OPTS                 <ARGS> Hand args over to the Java Virtual Machine.
                             In case of a heap-space error, assign more space with the
                             Arguments: "-Xms128g -Xmx128g" (allocates 128GB heap-space for the run)

.. _work:

Intermediate files
------------------

To run all processes separate from each other Nextflow creates intermediate files and directories
that are stored in the :file:`work` directory. You can delete the folder after the run::

    rm -fr work/

