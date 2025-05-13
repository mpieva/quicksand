.. _configuration-page:

Configuration
=============

.. _config:

Config-file
-----------

| **Create a custom config file**
| Most of the required flags won't change across quicksand processings, like the :code:`--db` or the :code:`--genomes` flag.
| To re-use parameters and other configs across runs, create a config-file e.g. :file:`quicksand.config` as described below and hand it over to the pipeline with the :code:`-c` flag:
::

    //
    // just put all the parameters that are used by default here

    params {
        db         = "path/to/kraken/Mito_db_kmer22"
        genomes    = "path/to/genomes"
        bedfiles   = "path/to/masked/"
    }

    //add singularity by default to the pipeline
    //bind: add shared drives (if necessary)
    //cacheDir: where to download the container

    singularity{
        enabled = true
        autoMounts = true
        runOptions = "--bind /mnt/"
        cacheDir   = "path/to/singularity/"
    }

Please see the available configuration options in the `nextflow documentation <https://www.nextflow.io/docs/latest/config.html#scope-singularity>`_

And run quicksand::

    nextflow run mpieva/quicksand -c nextflow.config [...] -profile singularity

Environmental variables
-----------------------

The following nextflow specific ENV variables can be set::

    NXF_SINGULARITY_CACHEDIR <path> Where to save the pulled Singularity-images
    NXF_WORK                 <path> Corresponds to the -w flag
    NXF_OPTS                 <ARGS> Hand args over to the Java Virtual Machine.
                             In case of a heap-space error, assign more space with the
                             Arguments: "-Xms10g -Xmx20g" (allocates 128GB heap-space for the run)

.. _work:

Intermediate files
------------------

Nextflow stores intermediate files and directories in the :file:`work` directory. You can delete the folder after the run::

    rm -fr work/

