.. _quickstart-page:
.. role:: bold

Quickstart
===========

| Use this section to test the :bold:`runability` of the pipeline
| A more comprehensive step-by-step-guide for the setup of the pipeline and the :bold:`test on real data` can be found in the :ref:`install-page` section.

Requirements
------------

See :ref:`requirements`. Make sure to have :code:`Nextflow` and :code:`Singularity or Docker` installed

Run quicksand
-------------

Test the runability of the quicksand pipeline using the included test-data::

	nextflow run mpieva/quicksand -profile singularity,test

| See the :code:`final_report.tsv` for a summary of the results. 
| See the :ref:`output` section for a detailed explaination of all the output files.               

Troubleshooting
---------------

- :bold:`heap-space error`: Nextflow related error, not enough memory allocated for the Java Virtual machine to run the processes. Manually set the memory for nextflow by typing :code:`export NXF_OPTS="-Xms128g -Xmx128g"` before the run
- :bold:`file-not found error`: Your Singularity installation might not allow auto-mounting of paths: See :ref:`singularity`  