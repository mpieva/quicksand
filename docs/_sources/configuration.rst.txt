.. _configuration-page:

Configuration
=============

.. _config:

Config-file
-----------

| **reduce the number of flags**
| Since most of the required flags won't change across runs - like the :code:`--db` or the :code:`--genome` flag, it is recommended to create a custom :file:`nextflow.config` file or a wrapper script to hand default values for these parameters over to the pipeline. For the config-file, write down all the static parameter/flags (see example below) and hand the file over to the pipeline with the :code:`-c` flag:

::

    nano ~/pipeline/nextflow.config

    --- GNU nano ---
    //add parameters to the pipeline
    params { 
        db         = "~/pipeline/data/kraken/Mito_db_kmer22"
        genome     = "~/pipeline/data/genomes"
        bedfiles   = "~/pipeline/data/masked/"
        report     = true
        analyze    = true
    }
    
    //add singularity-parameters to the pipeline
    //bind: add shared drives (if necessary)
    //cacheDir: where to download the container
    singularity{ 
        runOptions = "--bind /mnt/" 
        cacheDir   = "~/pipeline/singularity/"
    }

save the file by pressing :kbd:`ctrl` + :kbd:`x` and :kbd:`y`. Then run the pipeline with the following parameters::

    nextflow run mpieva/quicksand -c ~/pipeline/nextflow.config --split PATH/TO/INPUT/DATA

| **Global config-file**
| to avoid adding the :code:`-profile` flag to the command (in case one uses Docker or conda), it is possible to overwrite the default global config. Downloade the "docker config", add the desired local parameters and feed the file back to the pipeline with the :code:`-C` flag:

::

    nextflow config mpieva/quicksand -profile docker > ~/pipeline/nextflow.config
    nano ~/pipeline/nextflow.config

**Add** the parameters introduced above to the *params* scope. Then run the pipeline with::

    nextflow -C ~/piepline/nextflow.config  run mpieva/quicksand --split PATH/TO/SPLIT/DATA

.. attention::
    when overwriting the default configs, the :code:`-C nextflow.config` flag comes before the :code:`run mpieva/quicksand` command


Profiles
--------

| Profiles are a preinstalled set of parameters the pipeline starts with. 
| The following profiles are available:

    standard (no profile picked)   
        Run all processes of quicksand in a Singularity container. see the :ref:`singularity` section for more infos 
    conda      
        The non-containerized pipeline. All processes are executed with locally installed software and conda. Since the pipeline relies on some MPI-inhouse-packages that are not published on conda, it is not impossible to run the pipeline (see the :ref:`conda` section) with this profile outside the MPI EVA, but we highly recommend sticking to Singularity or Docker.
    docker
        Run all processes of quicksand in a Docker container. see the :ref:`singularity` section for more infos
    cluster       
        execute certain CPU-intensive processes on SGE cluster. **Be aware**: Because of the enourmous amount of intermediate files that nextflow creates and thus a lot of copying between shared drives, we recommend running quicksand without SGE on a local SSD-Drive (or even better: RAM-disc)!. Since experiences may vary, we still kept the option of running the pipeline on a cluster. 


Environmental variables
-----------------------

Environmental variables can be set to reduce the number of arguments handed over to the pipeline::

    Pipeline-variables
    QS_DB        <path>     Corresponds to the --db flag
    QS_GENOME    <path>     Corresponds to the --genome flag
    QS_BEDFILES  <path>     Corresponds to the --bedfiles flag
    QS_SPECMAP   <path>     Corresponds to the --bedfiles flag

    useful nextflow-variables
    NXF_SINGULARITY_CACHEDIR <path> Where to save the pulled Singularity-images
    NXF_WORK                 <path> Corresponds to the -w flag
    NXF_OPTS                 <ARGS> Hand args over to the Java Virtual Machine. 
                             In case of a heap-space error, assign more space with the
                             Arguments: "-Xms128g -Xmx128g" (allocates 128GB heap-space for the run)

Instead of working with :file:`nextflow.config` files, set the desired environmental variables::

    export QS_DB="~/pipeline/data/kraken/Mito_db_kmer22"
    export QS_GENOME="~/pipeline/data/genomes"
    export QS_BEDFILES="~/pipeline/data/masked/"
    export NXF_SINGULARITY_CACHEDIR="~/pipeline/singularity/"

And run the pipeline::
    
    nextflow run mpieva/quicksand --split PATH/TO/INPUT/DATA

.. _singularity:

Singularity
-----------

| To cite the Nextflow-docs [nxf_docs]_:
| "Singularity is a container engine alternative to Docker. The main advantages of Singularity is that it can be used with unprivileged permissions and doesn’t require a separate daemon process. These, along other features, like for example the support for autofs mounts, makes Singularity a container engine better suited the requirements of HPC workloads. Singularity is able to use existing Docker images, and pull from Docker registries."

So both the :code:`standard` and the :code:`docker` profiles use the same container to run the processes of the pipeline in. The image used to construct the containers is hosted on `dockerhub <https://hub.docker.com/repository/docker/merszym/quicksand>`_ and automatically pulled by the pipeline upon request (by running the pipeline). The Dockerfile used to create the image is hosted in the repository as well, see: :file:`quicksand/docker/Dockerfile`

.. attention::
    When running the pipeline with Singularity make sure that you:

        Bind the paths
            It happens that singularity is configured in a way that doesnt allow the automated mounting of paths into the container as intended by nextflow. If that is the case and files are either accessed or written to a shared drive, the process might crash and the pipeline exits with a **Path/File not found** error. If that is the case, make sure to add the path via the :code:`runOptions "--bind PATH"` to the *singularity* scope of your nextflow-config file in a way described in the :ref:`config section <config>`
        Set a cacheDir
            By default, nextflow downloads the Singularity image into the :file:`work/singularity` folder of your run. Since downloading and storing the image takes time and disc-space, it is recommended to set up a Singularity cacheDir. This can be done either by adding the :code:`cacheDir PATH` argument to the *singularity* scope of the custom :file:`nextflow.config` file (see :ref:`above <config>`). or by setting the environmental variable :code:`export NXF_SINGULARITY_CACHEDIR=PATH` before the run. With the directory set, Singularity will reuse your image.

.. [nxf_docs] https://www.nextflow.io/docs/latest/singularity.html


.. _conda:

Conda
-----

WIP


.. _work:

Intermediate files
------------------

To run all processes separate from each other Nextflow creates a lot of intermediate files and directories that are stored within the runDir in the :file:`work` directory. Since the file-tree inside the :file:`work` directory consumes much more disc-space than the whole "real" output of the pipeline, we recommend to delete the folder after the run. This can be either done over the command line::

    rm -fr work/

or via the :code:`nextflow clean` command::

    nextflow clean -f -q

to also remove hidden logs
