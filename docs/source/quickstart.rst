.. _quickstart-page:

Quickstart
===========

Requirements
------------

quicksand has two dependencies

:Nextflow: Version :code:`22.04` or above. `See here <https://www.nextflow.io/docs/latest/getstarted.html>`_
:Container: Please use `Singularity <https://sylabs.io/guides/3.0/user-guide/installation.html>`_ or `Docker <https://docs.docker.com/get-docker/>`_

.. tip::

    check the successful installation of the software by running::

        nextflow -v
            >>> nextflow version 22.10
        singularity --version
            >>> singularity version 3.7.2-dirty


Download the database
---------------------

THIS SECTION IS WORK IN PROGRESS

Download test-data
------------------

As input for the pipeline, download the Hominin "Hohlenstein-Stadel" mtDNA [1]_ into a directory `split`::

	wget -P split http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/BAM/mtDNA/HST.raw_data.ALL.bam


Run quicksand
-------------

nextflow pipelines can be executed directly from github. To run quicksand using the downloaded data-set type::

    nextflow run mpieva/quicksand -r v2.0 \
      -profile   singularity \
      --db       refseq/kraken/Mito_db_kmer22 \
      --genomes  refseq_rel220/genomes/ \
      --bedfiles refseq_rel220/masked/ \
      --split    split/


| The output of quicksand can be found in the directory **quicksand_v2.0/**
| See the :code:`final_report.tsv` for a summary of the results.
| See the :ref:`output-page` section for a detailed explaination of all the output files.

.. [1] http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/README