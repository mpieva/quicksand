Quickstart
===========

A comprehensive step-by-step-guide can be found in the :ref:`install-page` section.

Requirements
------------

Please ask your admin to install :code:`Nextflow` and :code:`Singularity`, otherwise see :ref:`requirements`

Create Datastructure
--------------------

Open your terminal with :kbd:`Win` + :kbd:`r` and type::
	
	mkdir quickstart && cd quickstart

	nextflow run mpieva/quicksand-build --outdir refseq

Run quicksand
-------------

For the run, we need some data: For that we download the Hohlenstein-Stadel mtDNA (please see the [README]_) ::
	
	wget -P split http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/BAM/mtDNA/HST.raw_data.ALL.bam

And run the quicksand pipeline::

	nextflow run mpieva/quicksand \
		--db 		refseq/kraken/Mito_db_kmer22 \
		--genome	refseq/genomes \
		--bedfiles	refseq/masked \
		--split 	split \
		--report 	\
		--analyze	\
		--byrg		

Please see the :ref:`output` section for an explaination of the output!

Problems?
---------
Some *classical* problems are:

- heap-space error: not enough heap-space allocated for the JVM. please type :code:`export NXF_OPTS="-Xms128g -Xmx128g"` before the run
- file-not found error: The Singularity installation doesnt allow auto-mounting of paths: See :ref:`singularity`  
- quicksand-build crashes: The indexing of the kraken-database requires a lot of RAM!


.. [README] http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/README


