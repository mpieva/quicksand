![MIT License](https://img.shields.io/github/license/mpieva/quicksand?style=for-the-badge)
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.11106450-ff69b4?style=for-the-badge)](https://zenodo.org/doi/10.5281/zenodo.11106450)

# quicksand

See [readthedocs](https://quicksand.readthedocs.io/en/latest/in_and_out.html) for the full documentation of the pipeline.

## Description

quicksand (**quick** analysis of **s**edimentary **an**cient **D**NA) is an open-source [Nextflow](https://doi.org/10.1038/nbt.3820) pipeline designed for rapid and accurate taxonomic classification of mammalian mitochondrial DNA (mtDNA) in aDNA samples. quicksand combines fast alignment-free classification using [KrakenUniq](https://doi.org/10.1186/s13059-018-1568-0) with downstream mapping ([BWA](https://github.com/mpieva/network-aware-bwa)), post-classification filtering, and ancient DNA authentication. quicksand is optimized for speed and portablity and requires either [Singularity](https://doi.org/10.1371/journal.pone.0177459) or [Docker](https://www.docker.com/).

## Workflow

<p align=center>
    <img src="assets/docs/workflow_v2.3.png" alt="Graphical representation of the pipeline workflow" width='800px'>
</p>

## Quickstart

### Requirements

To run Nextflow, you need a POSIX-compatible system (e.g., Linux or macOS). quicksand was developed and tested on Linux (x86_64 architecture)

To run quicksand, please install

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) v22.10 or larger
- [Singularity](https://sylabs.io/singularity/) or [Docker](https://www.docker.com/)

**Note:** To run quicksand in singularity, your kernel needs to support user-namespaces (see [here](https://github.com/apptainer/singularity/issues/5240#issuecomment-618405898) or [here](https://github.com/apptainer/singularity/issues/6341)).


### Prepare Input

The input for quicksand is a directory with user-supplied files in BAM or FASTQ format. Adapter-trimming, overlap-merging and sequence demultiplexing need to be performed by the user prior to running quicksand. Provide the directory with the `--split` flag

> [!CAUTION]
> Each input-file should correspond to a single sequence-library. The processing of merged libraries with quicksand can lead to sequence loss because of the PCR-deduplication step with bam-rmdup 

#### Download Test-file

As a test file, download the Hohlenstein-Stadel mtDNA (please see the [README](http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/README) for more information)

```bash
wget -P split \
http://ftp.eva.mpg.de/neandertal/Hohlenstein-Stadel/BAM/mtDNA/HST.raw_data.ALL.bam
```

### Create Reference Database

The required KrakenUniq database, the reference genomes for mapping and the bed-files for low-complexity filtering are available on the MPI EVA FTP Servers. Custom versions of the reference material can be created with the [quicksand-build pipeline](https://github.com/mpieva/quicksand-build)

#### Create Test Database

For the quickstart of quicksand, create a fresh database containing only the Hominidae mtDNA reference genomes (runtime: ~3-5 minutes)

```bash
nextflow run mpieva/quicksand-build -r v3.1 \
  --include  Hominidae \
  --outdir   refseq \
  -profile   singularity
```

#### Download Full Database

 To download the full reference database (~60GB), use this command:

```bash
latest=$(curl http://ftp.eva.mpg.de/quicksand/LATEST)
wget -r -np -nc -nH --cut-dirs=3 --reject="*index.html*" -q --show-progress -P refseq http://ftp.eva.mpg.de/quicksand/build/$latest
```

### Run quicksand

quicksand is executed directly from github. With the databases created and the testdata downloaded, run the pipeline as follows:

```bash
# set this if you encounter a heap-space error to increase the memory that is used by nextflow
export NXF_OPTS="-Xms10g -Xmx15g" # increase or decrease the numbers as required

nextflow run mpieva/quicksand -r v2.4 \
  --db        refseq/kraken/Mito_db_kmer22/ \
  --genomes   refseq/genomes/ \
  --bedfiles  refseq/masked/ \
  --split     split/ \
  -profile    singularity #mind the single dash!
```

### Output

Please see the [documentation](https://quicksand.readthedocs.io/en/latest/in_and_out.html) for a comprehensive description of the output!

### Common Errors
A collection of common nextflow-errors and how to solve them

#### Heap Space
```
 -- Check '.nextflow.log' file for details
ERROR ~ Java heap space

 -- Check '.nextflow.log' file for details
ERROR ~ Execution aborted due to an unexpected error
```

Heap space errors can occur if nextflow itself requires more memory than provided by default (e.g. when screening too many samples in parallel). You can increase the heap-space as needed (e.g., to 5gb) with 

```
export NXF_OPTS="-Xms5g -Xmx5g"
``` 

## References and Citation

If you use quicksand in your research, please cite the quicksand-preprint as follows:

> Szymanski, Merlin, Johann Visagie, Frederic Romagne, Matthias Meyer, and Janet Kelso. 2025.
> “Quick Analysis of Sedimentary Ancient DNA Using Quicksand”,
> [https://doi.org/10.1101/2025.08.01.668088](https://doi.org/10.1101/2025.08.01.668088).

This pipeline uses code inspired by the [nf-core](https://nf-co.re) initative, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> Ewels, Philip A., Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso, and Sven Nahnsen. 2020.
> “The Nf-core Framework for Community-curated Bioinformatics Pipelines”.
> Nature Biotechnology 38 (3): 276–78. [https://doi.org/10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x).
