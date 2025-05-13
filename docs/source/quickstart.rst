.. _quickstart-page:

Quickstart
===========

Requirements
------------

quicksand has two dependencies

:Nextflow: Version :code:`22.04` or above. `See here <https://www.nextflow.io/docs/latest/getstarted.html>`_
:Containerization-Software: Please use `Singularity <https://sylabs.io/guides/3.0/user-guide/installation.html>`_ or `Docker <https://docs.docker.com/get-docker/>`_

.. tip::

    check the successful installation of the software by running::

        nextflow -v
            >>> nextflow version 22.10
        singularity --version
            >>> singularity version 3.7.2-dirty


Download test-data
------------------

The input for quicksand is a directory with user-supplied files in BAM or FASTQ format. 
Adapter-trimming, overlap-merging and sequence demultiplexing need to be performed by the user prior to running quicksand. 
Provide the directory with the :code:`--split DIR` flag. 

As input for the quickstart, download the Hominin "Hohlenstein-Stadel" mtDNA [1]_ into a directory `split`::

	wget -q --show-progress -P split http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/BAM/mtDNA/HST.raw_data.ALL.bam


Download the database
---------------------

The required KrakenUniq database, the reference genomes for mapping and the bed-files for low-complexity filtering are available on the 
MPI EVA FTP Servers. Custom versions of the reference material can be created with the quicksand-build pipeline

For quickstarting quicksand, create a fresh database containing only the Hominidae mtDNA reference genomes (runtime: ~3-5 minutes)::

nextflow run mpieva/quicksand-build -r v3.0 \
  --include  Hominidae \
  --outdir   refseq \
  -profile   singularity


Alternatively, download the most full datastructure from the MPI EVA FTP SERVERS (~50 GB)::

    latest=$(curl http://ftp.eva.mpg.de/quicksand/LATEST)
    wget -r -np -nc -nH --cut-dirs=3 --reject="*index.html*" -q --show-progress -P refseq http://ftp.eva.mpg.de/quicksand/build/$latest


Run quicksand
-------------

quicksand is executed directly from github. With the databases created and the testdata downloaded, run the pipeline as follows::

    nextflow run mpieva/quicksand -r v2.4 \
      -profile   singularity \
      --db       refseq/kraken/Mito_db_kmer22 \
      --genomes  refseq/genomes/ \
      --bedfiles refseq/masked/ \
      --split    split/


The output of quicksand can be found in the directory **quicksand_v2.4/**

See the :code:`final_report.tsv` and :code:`filtered_report_0.5p_0.5b.tsv` for a summary of the results.

See the :ref:`output-page` section for a detailed explaination of all the output files.

.. [1] http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/README