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

Run the pipeline::

    nextflow run mpieva/quicksand -c nextflow.config --split path/to/split -profile singularity



Environmental variables
-----------------------

The following nextflow specific ENV variables can be set::

    NXF_SINGULARITY_CACHEDIR <path> Where to save the pulled Singularity-images
    NXF_WORK                 <path> Corresponds to the -w flag
    NXF_OPTS                 <ARGS> Hand args over to the Java Virtual Machine. 
                             In case of a heap-space error, assign more space with the
                             Arguments: "-Xms128g -Xmx128g" (allocates 128GB heap-space for the run)

.. _singularity:

Singularity
-----------

| To cite the `Nextflow docs <https://www.nextflow.io/docs/latest/singularity.html>`_:
| "Singularity is a container engine alternative to Docker. The main advantages of Singularity is that it can be used with 
| unprivileged permissions and doesn’t require a separate daemon process. These, along other features, like for example 
| the support for autofs mounts, makes Singularity a container engine better suited the requirements of HPC workloads. 
| Singularity is able to use existing Docker images, and pull from Docker registries."

The containers used in the pipeline are explained here: :ref:`container` 

.. attention::
    When running the pipeline with Singularity make sure that you:

        Bind the paths
            It happens that singularity is configured in a way that doesnt allow the automated mounting of paths into the container as 
            intended by nextflow. If that is the case and files are either accessed or written to a shared drive, 
            the process might crash and the pipeline exits with a **Path/File not found** error. If that is the case, make sure to 
            add the path via the :code:`runOptions "--bind PATH"` to the *singularity* scope of your nextflow-config file in a way 
            described in the :ref:`config section <config>`
        Set a cacheDir
            By default, nextflow downloads the Singularity image into the :file:`work/singularity` folder of your run. 
            Since downloading and storing the image takes time and disc-space, it is recommended to set up a Singularity cacheDir. 
            This can be done either by adding the :code:`cacheDir PATH` argument to the *singularity* scope of the custom :file:`nextflow.config` file 
            (see :ref:`above <config>`). or by setting the environmental variable :code:`export NXF_SINGULARITY_CACHEDIR=PATH` before the run. 
            With the directory set, Singularity will reuse the pulled images.


.. _work:

Intermediate files
------------------

To run all processes separate from each other Nextflow creates intermediate files and directories 
that are stored in the :file:`work` directory. Since the file-tree inside the :file:`work` directory consumes much disc-space, 
we recommend to delete the folder after the run::

    rm -fr work/

