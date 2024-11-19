![MIT License](https://img.shields.io/github/license/mpieva/quicksand?style=for-the-badge)
![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.11106450-ff69b4?style=for-the-badge&link=https://zenodo.org/doi/10.5281/zenodo.11106450)

<h1 style="border:0px;padding-bottom:0px;margin-bottom:0px">quicksand</h1>
<p style="color:grey;border-bottom:1px solid lightgrey">quick analysis of sedimentary ancient DNA</p>

See the [documentation](https://quicksand.readthedocs.io/en/latest/in_and_out.html) for a comprehensive documentation of the pipeline.

## Description

quicksand is a bioinformatic pipeline for the analysis and taxonomic binning of (target enriched) ancient, mitochondrial, sedimentary DNA. quicksand uses [krakenuniq](https://doi.org/10.1186/s13059-018-1568-0) for metagenomic classification, [BWA](https://github.com/mpieva/network-aware-bwa) for the mapping of DNA sequences and analyses mapped sequences for DNA deamination patterns.

Optimized for speed and portablity, quicksand is written in [Nextflow](https://doi.org/10.1038/nbt.3820) and requires either [Singularity](https://doi.org/10.1371/journal.pone.0177459) or [Docker](https://www.docker.com/).

## Workflow

<p align=center>
    <img src="assets/docs/v1.2.png" alt="Graphical representation of the pipeline workflow" width='800px'>
</p>

## Quickstart

### Requirements

To run the pipeline, please install

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) v22.10 or larger
- [Singularity](https://sylabs.io/singularity/) or [Docker](https://www.docker.com/)

### Input

The pipeline accepts demultiplexed, adapter-trimmed and overlap-merged bam and fastq files. Put all files in one directory, name the files `DIR/{READGROUP}.{bam, fastq}`. Provide the directory with the `--split` flag

### Download Datastructure

To run quicksand a kraken database for metagenomics classification, the reference genomes for mapping and a set of bed-files are required for the run of the pipeline.

For the most recent RefSeq releases please download the quicksand-datastructure
here:

```bash
latest=$(curl http://ftp.eva.mpg.de/quicksand/LATEST)
wget -r -np -nc -nH --cut-dirs=3 --reject="*index.html*" -q --show-progress -P refseq http://ftp.eva.mpg.de/quicksand/build/$latest
```

This step takes a while! Make yourself a coffee and relax

For a custom creation of the datastructure see the [quicksand-build pipeline](https://github.com/mpieva/quicksand-build)

### Download Test-data

To run quicksand with real data, download the Hohlenstein-Stadel mtDNA (please see the [README](http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/README) for more information) as input

```bash
wget -P split \
http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/BAM/mtDNA/HST.raw_data.ALL.bam
```

### Run quicksand

quicksand is executed directly from github, no local build is required. With the databases and the testdata downloaded, run the pipeline.

```bash
nextflow run mpieva/quicksand -r v2.3 \
  --db        refseq/kraken/Mito_db_kmer22/ \
  --genomes   refseq/genomes/ \
  --bedfiles  refseq/masked/ \
  --split     split/ \
  -profile    singularity
```

### Output

Please see the [documentation](https://quicksand.readthedocs.io/en/latest/in_and_out.html) for a comprehensive description of the output!

## References

This pipeline uses code inspired by the [nf-core](https://nf-co.re) initative, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.
