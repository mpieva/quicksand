# Quicksand

![MIT License](https://img.shields.io/github/license/mpieva/quicksand?style=for-the-badge)
![Status Warning](https://img.shields.io/badge/status-WIP-ff69b4?style=for-the-badge)

**Quicksand**: quick analysis of sedimentary ancient DNA \
Please see https://mpieva.github.io/quicksand for a **comprehensive documentation**.

## Description
quicksand is a bioinformatic pipeline for the analysis and taxonomic binning of target enriched ancient mammalian mitochondrial sedimentary DNA.
With the general workflow and background described in [_Slon et al., 2017_](https://science.sciencemag.org/content/sci/suppl/2017/04/26/science.aam9695.DC1/aam9695_SM.pdf), 
quicksand uses [krakenuniq](https://doi.org/10.1186/s13059-018-1568-0) for the metagenomic classification of reads and is optimized for speed and reproducibility. 
The pipeline is written in [Nextflow](https://doi.org/10.1038/nbt.3820) and requires [Singularity](https://doi.org/10.1371/journal.pone.0177459) or [Docker](https://www.docker.com/).

While the default settings are optimized and tested for the assignment of mammalian mtDNA, `quicksand` can be used together with databases constructed from the whole RefSeq mtDNA database
(see [HERE](https://www.github.com/mpieva/quicksand-build) for more information about the database construction).

#### Workflow   

<p align=center>
    <img src="assets/pipeline_overview.png" alt="Graphical overview over the pipeline workflow" width=500 height=500>
</p>

#### Input

The pipeline starts with genomic data (bam and fastq-files)
Put all your `READGROUP.bam` and all your `READGROUP.{fq,fq.gz,fastq,fastq.gz}` files into one directory and use this directory as input
with the `--split` statement 

**Note #1:**
- In a perfect world, the files are all _demultiplexed_, _adapter-trimmed_ and _overlap-merged_
- the pipeline interprets paired-end fastq files as single-read files upon conversion to bam format (so merge them, or accept it)
- paired-end BAM files are filtered out by default 

**Note #2:**
The pipeline includes a splitBam process, however, that one will not work with libraries and index-combinations that were not produced by the [MPI EVA CoreUnit](https://www.eva.mpg.de/genetics/index/).
So make sure to demultiplex your input files before processing (as it is usually done by the Illumina sequencers) 

#### Output

- For each readgroup: The processed and binned reads in `.bam`-format at each stage of the pipeline (see workflow above). 
- Summary stats for each readgroup and assigned family: Number of assigned reads, mapped reads, unique reads, bedfiltered reads, deaminated reads and some more basic stats 

## Quickstart
### Requirements

To run the pipeline, please install
 - Nextflow
 - Singularity

### Create Datastructure

But WAIT! to make a metagenomic classification, we need a reference database, some reference genomes and a taxonomy...
To run the `quicksand` pipeline, we have to run the supplementary pipeline `quicksand-build` in advance to do
exactly that for us. This pipeline will download the taxonomy from NCBI/taxonomy, the mitochondrial genomes from NCBI/RefSeq
and build the kraken-database with the specified settings. 

For this test we will create a database containing only the _Primates_. If you want to use the pipeline for the analysis
of all mammalian or _everything_ in RefSeq, please see the [README](https://www.github.com/mpieva/quicksand-build) of this pipeline

To run the pipline open your terminal and type:

```bash
    mkdir quickstart && cd quickstart
    nextflow run mpieva/quicksand-build --outdir refseq --include Primates
```

And wait. Especially the download of the taxonomy takes quite long
    
### Run quicksand

With the databases created in `refseq` we can now run the actual pipeline.
For the run itself, we need some data: So we download the Hohlenstein-Stadel mtDNA (please see the [README](http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/README) for more information)

```bash
    wget -P split http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/BAM/mtDNA/HST.raw_data.ALL.bam
```

Then we can finally run the pipeline:

```bash
    nextflow run mpieva/quicksand \
        --db        refseq/kraken/Mito_db_kmer22 \
        --genomes   refseq/genomes \
        --bedfiles  refseq/masked \
        --split     split \
        -profile    singularity
```

After running the pipeline, please see the `final_report.tsv` for a summary of the results.
Please see the [DOCS](https://mpieva.github.io/quicksand/usage.html#output) for a comprehensive description of the output!

In case of questions or bugs, 
feel free to contact me or open an issue

