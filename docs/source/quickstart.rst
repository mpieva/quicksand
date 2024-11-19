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

Download the most recent datastructure for running the quicksand pipeline here::

    latest=$(curl http://ftp.eva.mpg.de/quicksand/LATEST)
    wget -r -np -nc -nH --cut-dirs=3 --reject="*index.html*" -q --show-progress -P refseq http://ftp.eva.mpg.de/quicksand/build/$latest

This step takes a while! Make yourself a coffee and relax

Download test-data
------------------

As input for the pipeline, download the Hominin "Hohlenstein-Stadel" mtDNA [1]_ into a directory `split`::

	wget -q --show-progress -P split http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/BAM/mtDNA/HST.raw_data.ALL.bam


Run quicksand
-------------

nextflow pipelines can be executed directly from github. To run quicksand using the downloaded data-set type::

    nextflow run mpieva/quicksand -r v2.3 \
      -profile   singularity \
      --db       refseq/kraken/Mito_db_kmer22 \
      --genomes  refseq/genomes/ \
      --bedfiles refseq/masked/ \
      --split    split/


| The output of quicksand can be found in the directory **quicksand_v2.3/**
| See the :code:`final_report.tsv` and :code:`filtered_report_0.5p_0.5b.tsv` for a summary of the results.
| See the :ref:`output-page` section for a detailed explaination of all the output files.

.. [1] http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/README