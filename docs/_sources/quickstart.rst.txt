.. role:: bold
.. _quickstart-page:

Quickstart
===========

This quickstart section is slightly more detailed than 
the one in the `README <https://github.com/mpieva/quicksand/blob/master/README.md>`_. of the quicksand repository. A more comprehensive 
step-by-step-guide for the setup can be found in the :ref:`install-page` section.

Requirements
------------

Please make sure you have installed 

* :code:`Nextflow` and 
* :code:`Singularity` or :code:`Docker` (see :ref:`requirements`) 

Create datastructure
--------------------

The required underlying datastructure of the pipeline is described in the :ref:`quicksand_build-page` section

In short: You need a :bold:`reference database`, the :bold:`reference genomes` and a :bold:`taxonomy`. 
Use the supplementary pipeline :code:`quicksand-build` (once) to download the taxonomy from NCBI/taxonomy, all mitochondrial
genomes from NCBI/RefSeq and create the required databases and files.

For this quickstart-session create a datastructure containing only the :bold:`Primates`::
	
	mkdir quickstart && cd quickstart
	nextflow run mpieva/quicksand-build --outdir refseq --include Primates

| :bold:`Note:`
| Building the database requires ~40G of RAM
| Be patient, downloading the taxonomy plus the creation of the database might take :bold:`~1h`.

Run quicksand
-------------

| With the datastructure created in the :code:`refseq` directory, the pipeline is ready to be used.
| As :bold:`input` for the pipeline, download the Hominin "Hohlenstein-Stadel" mtDNA [1]_ into a directory :bold:`split` 
::
	
	wget -P split http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/BAM/mtDNA/HST.raw_data.ALL.bam

And run the quicksand pipeline::

	nextflow run mpieva/quicksand \
	   --db        refseq/kraken/Mito_db_kmer22 \
	   --genome    refseq/genomes \
	   --bedfiles  refseq/masked \
	   --split     split \
	   -profile    singularity

| See the :code:`final_report.tsv` for a summary of the results. 
| See the :ref:`output` section for a detailed explaination of all the output files.               


Filter the Results
------------------

As can be seen in the :code:`final_report.tsv`, not all sequences were assigned to Homindae, but to a couple of other Primate families as well.
The assignment of false positive taxa is a well-known problem of kmer-based assignment methods and additional filters need to be applied.

Based on simulated data, our recommended cutoffs are:

- :bold:`FamPercentage` cutoff of 1% and/or
- :bold:`ProportionMapped` cutoff of 0.5-0.7.

The kmer-information is also indicative. If the :bold:`FamilyKmers` and :bold:`KmerCoverage` values are low and
the :bold:`KmerDupRate` value is high, the assigment of the family is only based on a small number of kmers within the reads 

Troubleshooting
---------------

- :bold:`heap-space error`: Nextflow related error, not enough memory allocated for the Java Virtual machine to run the processes. Manually set the memory for nextflow by typing :code:`export NXF_OPTS="-Xms128g -Xmx128g"` before the run
- :bold:`file-not found error`: Your Singularity installation might not allow auto-mounting of paths: See :ref:`singularity`  

.. [1] http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/README