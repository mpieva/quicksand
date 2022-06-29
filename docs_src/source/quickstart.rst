Quickstart
===========

This quickstart section can also be found in the README of the quicksand repository.
A more comprehensive step-by-step-guide can be found in the :ref:`install-page` section.

Requirements
------------

Please ask your admin to install :code:`Nextflow` and :code:`Singularity`, otherwise see :ref:`requirements`

Create Datastructure
--------------------

To make a metagenomic classification, we need a reference database, reference genomes and a taxonomy.
To run the :code:`quicksand` pipeline, we have to run the supplementary pipeline :code:`quicksand-build` (once) in advance to do
exactly that for us. The :code:`quicksand-build` pipeline will download the taxonomy from NCBI/taxonomy, the mitochondrial genomes from NCBI/RefSeq
and build the kraken-database with the specified settings.

For the quickstart-session we create a database containing only the _Primates_ by providing the :code:`--include Primates` flag. If you 
want to use the pipeline for the analysis of the whole mammalian diversity use :code:`--include Mammalia` or remove the flag if 
the database should contain _everything_ in RefSeq. 

To run the pipline open your terminal and type:

Open your terminal with :kbd:`Win` + :kbd:`r` and type::
	
	mkdir quickstart && cd quickstart
	nextflow run mpieva/quicksand-build --outdir refseq --include Primates

Please be patient, the download of the taxonomy and the creation of the database might take ~1h
After that, we are ready to run the :code:`quicksand` pipeline

Run quicksand
-------------

With the databases created in :code:`refseq` we can now run the actual pipeline.
However, before that we need some data: So download the Hohlenstein-Stadel mtDNA (please see the [README]_ for more information) ::
	
	wget -P split http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/BAM/mtDNA/HST.raw_data.ALL.bam

And run the quicksand pipeline::

	nextflow run mpieva/quicksand \
	   --db        refseq/kraken/Mito_db_kmer22 \
	   --genome    refseq/genomes \
	   --bedfiles  refseq/masked \
	   --split     split \
           -profile    singularity \
           --byrg


After running the pipeline, please see the :code:`final_report.tsv` for a summary of the results.                
Please see the :ref:`output` section for an explaination of the output!

Problems?
---------
Some *classical* problems are:

- heap-space error: not enough heap-space allocated for the JVM. please type :code:`export NXF_OPTS="-Xms128g -Xmx128g"` before the run
- file-not found error: The Singularity installation doesnt allow auto-mounting of paths: See :ref:`singularity`  
- quicksand-build crashes: The indexing of the kraken-database requires a lot of RAM!


.. [README] http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/README


