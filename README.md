<h1 style="border:0px;padding-bottom:0px;margin-bottom:0px">Quicksand</h1>
<p style="color:grey;border-bottom:1px solid lightgrey">A quick analysis of sedimentary ancient DNA</p>

![Singularity](https://img.shields.io/badge/run_with-Singularity-ff69b4?style=for-the-badge)
![Docker](https://img.shields.io/badge/run_with-Docker-0db7ed?style=for-the-badge)
![MIT License](https://img.shields.io/github/license/mpieva/quicksand?style=for-the-badge)


See the [Github Pages](https://mpieva.github.io/quicksand) for a comprehensive documentation of the pipeline.

<!-- TOC -->
- [Description](#description)
    - [Workflow](#workflow)
    - [Input](#input)
    - [Output](#output)
- [Quickstart](#quickstart)
    - [Requirements](#requirements)
    - [Create Datastructure](#create-datastructure)
    - [Run quicksand](#run-quicksand)
<!-- /TOC -->

## Description
quicksand is a bioinformatic pipeline for the analysis and taxonomic binning of (target enriched) ancient, mitochondrial, sedimentary DNA.

With the workflow and background described in [_Slon et al., 2017_](https://science.sciencemag.org/content/sci/suppl/2017/04/26/science.aam9695.DC1/aam9695_SM.pdf), quicksand uses [krakenuniq](https://doi.org/10.1186/s13059-018-1568-0) and [BWA](https://github.com/mpieva/network-aware-bwa) for the metagenomic classification and mapping of reads. Quicksand is optimized for speed and reproducibility. The pipeline is written in [Nextflow](https://doi.org/10.1038/nbt.3820) and requires either [Singularity](https://doi.org/10.1371/journal.pone.0177459) or [Docker](https://www.docker.com/).

While the default settings are optimized and tested for the assignment of mammalian mtDNA, quicksand can be combined with databases constructed from the whole RefSeq mtDNA database
(see [HERE](https://www.github.com/mpieva/quicksand-build)).

### Workflow   

<p align=center>
    <img src="assets/pipeline_overview_v1.5.png" alt="Graphical overview over the pipeline workflow" width='800px'>
</p>

### Input

The pipeline accepts `bam` and `fastq` files.\
Collect all files named `READGROUP.bam` and/or `READGROUP.{fq,fq.gz,fastq,fastq.gz}` into one directory.  

**Notes**
- The files should be _demultiplexed_, _adapter-trimmed_ and _overlap-merged_
- `fastq` files are converted to single-read `bam` files (in case the input is paired-end)
- paired-end reads are filtered from the `bam` files by default 

**Note #2:**\
The pipeline includes a splitBam process for demultiplexing, however, it is restricted to libraries and index-combinations produced by the [MPI EVA CoreUnit](https://www.eva.mpg.de/genetics/index/).

### Output
- For each readgroup: quicksand outputs the processed and binned reads in `.bam`-format at each stage of the pipeline (see workflow above). 
- Summary stats for each readgroup and assigned family: Number of assigned reads, mapped reads, unique reads, bedfiltered reads, deaminated reads.

## Quickstart
### Requirements

To run the pipeline, please install
 - Nextflow
 - Singularity or Docker

And create the underlying datastructure

### Create Datastructure

To make a metagenomic classification, a reference database, some reference genomes and the taxonomy is required.

To run quicksand, execute the supplementary pipeline `quicksand-build`  in advance to do exactly that. This pipeline will download the taxonomy from NCBI/taxonomy, the mitochondrial genomes from NCBI/RefSeq
and build the kraken-database with the specified settings. 

For this README, create a database containing only the _Primate mtGenomes_. To use the pipeline for the analysis of all mammalian mtGenomes or _everything_ in RefSeq, please see the [README](https://www.github.com/mpieva/quicksand-build) of this pipeline

To run the pipline open your terminal and type:

```bash
    mkdir quickstart && cd quickstart
    nextflow run mpieva/quicksand-build --outdir refseq --include Primates
```

And wait. Especially the download of the taxonomy takes ~1h
    
### Run quicksand

With the databases created in `refseq` we can now run the actual pipeline.
To do that, download the Hohlenstein-Stadel mtDNA (please see the [README](http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/README) for more information) as input

```bash
    wget -P split http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/BAM/mtDNA/HST.raw_data.ALL.bam
```

Then run the pipeline:

```bash
    nextflow run mpieva/quicksand \
        --db        refseq/kraken/Mito_db_kmer22 \
        --genomes   refseq/genomes \
        --bedfiles  refseq/masked \
        --split     split \
        -profile    singularity
```

After running the pipeline, please see the `final_report.tsv` for a summary of the results. Please see the [docs](https://mpieva.github.io/quicksand/usage.html#output) for a comprehensive description of the output!

In case of questions or bugs, \
feel free to contact me or open an issue
