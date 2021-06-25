Description
===========

This repository contains an mammlian mitochondrial ancientDNA analysis pipeline. 
With the overall workflow of the pipeline described in Slon et al., 2017, this 
implementation of the pipeline uses KrakenUniq for classification and modifies several steps. 
The pipeline is written in Nextflow and makes use of Singularity for reproducibility.

The Input of the pipeline consists of demultiplexed, trimmed and overlap-merged reads in either .bam or .fastq format. 
For each readgroup, the pipeline outputs the raw reads in .bam-format at each stage of the pipeline. 
Additionally the number of reads extracted, mapped, deduplicated, and bedfiltered as well as the 
estimated DNA-damage is reported for a quick overview in one big summary-file. See the Test 
section for an overview over the output-files. Taxonomic assignments are evaluated on a family level
