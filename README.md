# Sediment pipeline (Nextflow edition)

A pipeline to preprocess sequencing data and output classified reads mapped against the best fitting reference genome

**input:** The Pipeline includes an optional MPI-specific demultiplexer. The accepted input is thus either a bam/indexfile pair or already demultiplexed files (both bam, fastq files allowed).
**output:** The pipeline ends with classified reads mapped against the corresponding reference genome

This pipeline is a reimplmentation of the MPI Sediment processing pipeline described in Slon et al., 2017 using Kraken for classification, 
and to implement some other requested changes. The pipeline is based on Nextflow to provide a reproducible and distributable pipeline.

This README is split into two sections. The first one covers the use of the pipeline **inside** the Max-Planck-Institute with
access to the inhouse scripts and packages. The second part covers the more complex setup outside the institute including Docker

# Section 1 - MPI users
## Prerequisites

You need to have nextflow and conda installed together. If `nextflow` is not available by default, create (if not done before) and activate an environment that contains both Nextflow and the `conda` tool itself:

```bash
conda create -n nextflow -c bioconda nextflow conda
conda activate nextflow
```

You'll see that `(nextflow)` appears in front of your shell prompt when the environment is active. The environment will stay active until you close the shell session or use `conda deactivate`. The next time you log in, you will have to activate the environment again using `conda activate nextflow`.


## Quickstart

Given the input files BAM,RG or the SPLIT-DIR, the easiest way to run the pipeline is to use a wrapper-script.
Download the script:
```bash
git clone https://vcs.eva.mpg.de/merlin_szymanski/nf_wrap
```
Execute the wrapper-script with: 
```bash
bash nf_wrap/run.sh --bam BAM --rg RG
``` 
or
```bash
bash nf_wrap/run.sh --split DIR
``` 

## Manual Setup (Working around two bugs in Nexflow 😔)


If you have **not used the wrapper before** (described above), you'll have to create a small configuration file in your home directory before you're able to run a Nextflow pipeline straight from this repository.

### Preparing to run the pipeline
Activate the environment (see the prerequisites)

```bash
conda activate nextflow
```

Next, you'll need to create a configuration file `~/.nextflow/scm` that configures `vcs.eva` as a source for Nextflow pipelines. You can do so by cutting and pasting the following line into your shell. _You only have to do this once._
 
```bash
mkdir -p "${HOME}/.nextflow" && curl -s https://vcs.eva.mpg.de/visagie/sediment_nf/snippets/8/raw >"${HOME}/.nextflow/scm"
```

You can now use `nextflow pull` to pull this workflow from `vcs.eva` into your account. It is written into a special cache directory inside your home directory, so you do not have to be in any special location in the filesystem when you do this:
   
```bash
nextflow pull visagie/sediment_nf -hub eva
```

When the pipeline is updated in future, Nextflow will notify you that a newer version is available. At that point, you can update to the most recent version with the following command:
   
```bash
nextflow pull visagie/sediment_nf
```

### Running the pipeline

The pipeline requires a number of flags to be present. See the **Flags** -section for further documentation

```bash
nextflow run visagie/sediment_nf --db </path/to/kraken.db> { --bam <input bamfile> --rg <index file> | --split <path/to/split-dir> } --genome <path/to/reference/database> --bedfiles <path/to/bedfiles> [--dedup, --keeppaired, --filterunmapped, --specmap <specmap-file>, ]
```

If you wish to *resume* a pipeline run (e.g. if you stopped it for some reason, and you do not want it to redo the steps it had already completed), you need to add the flag `-resume` (note: just one `-`) after `run`:

```bash
nextflow run -resume visagie/sediment_nf --db </path/to/kraken.db> { --bam <input bamfile> --rg <index file> | --split <path/to/split-dir> } --genome <path/to/reference/database> --bedfiles <path/to/bedfiles> [--dedup, --keeppaired, --filterunmapped, --specmap <specmap-file>, ]
```

# Section 2 - Outside the Institute

To run the pipeline outside the institute docker is required to set up an appropriate environment
you can inspect the image at:

https://hub.docker.com/r/merszym/sediment_nf

## Prerequisites

1. Nextflow
2. Docker

## Quickstart/Test

If you have `nextflow` and `docker` installed, you can **test** the pipeline
clone this repository:
```bash
git clone https://vcs.eva.mpg.de/visagie/sediment_nf
```
activate (or create) an environment that contains nextflow:
```bash
conda create -n nextflow -c bioconda nextflow
conda activate nextflow
```
Unpack the Kraken-Test-Database (too large for github) and run the test-set in a separate folder
```bash
mkdir test_run && cd test_run
tar -xvzf../sediment_nf/assets/test/kraken/database.tar.gz
nextflow run ../sediment_nf/main.nf --split ../sediment_nf/assets/test/split/ --genome ../sediment_nf/assets/test/genomes/ --bedfiles ../sediment_nf/assets/test/masked/ --db TestDB/ --specmap ../sediment_nf/assets/test/genomes/specmap.tsv -profile docker
```

Be careful, the `-profile docker` has only one dash!

## Set up the pipeline
 
With the test-run being successful the next step is to set up the required Databases.
That process is described here:

https://github.com/MerlinSzymanski/datastructure_nf/

## Running the pipeline
(TODO: where is the pipeline hosted?)

The pipeline requires a number of flags to be present. See the **Flags** -section for further documentation

```bash
nextflow run sediment_nf/main.nf --db </path/to/kraken.db> { --bam <input bamfile> --rg <index file> | --split <path/to/split-dir> } --genome <path/to/reference/database> --bedfiles <path/to/bedfiles> [--dedup, --keeppaired, --filterunmapped, --specmap <specmap-file>, ]
```

If you wish to *resume* a pipeline run (e.g. if you stopped it for some reason, and you do not want it to redo the steps it had already completed), you need to add the flag `-resume` (note: just one `-`) after `run`:

```bash
nextflow run -resume sediment_nf/main.nf --db </path/to/kraken.db> { --bam <input bamfile> --rg <index file> | --split <path/to/split-dir> } --genome <path/to/reference/database> --bedfiles <path/to/bedfiles> [--dedup, --keeppaired, --filterunmapped, --specmap <specmap-file>, ]
```

# Flags

### required arguments:
**Either**
Flag | Input Type | Description
--- | --- | ---
--bam | FILE | Still multiplexed BAM file
--rg | FILE | Tab-separated file containing index combinations. Format: 'LibID<tab>P7<tab>P5'
  
**Or**
Flag | Input Type | Description
--- | --- | ---
--split | DIR | Directory with already split bamfiles or demultiplexed, merged fastq-files. files should be named: 'Readgroup.(bam,fastq)' 

**And**
Flag | Input Type | Description
--- | --- | ---
--db | DIR | Path to the Kraken database to use
--genome | DIR | Path to reference genomes for BWA-Mappings
--bedfiles | DIR | Path to bedfiles masking the genomes for
      
### optional arguments:
Flag | Input Type | Description
--- | --- | ---
--specmap | FILE | Config-file. Force BWA-mappings of a family to assigned species, overwrite Krakens species-level assignment. Format: 'Family<tab>Species_name,Species_name'. Species_name must correspond to the filename in the genomes diretory
--cutoff | N | Length cutoff after BWA-Mapping (default: 35)
--quality | N | Mapping Quality filter after BWA Mapping (default: 25)
--bwacutoff | N | Cutoff for number of kraken matches (default: 0)
--keeppaired | - | Keep paired reads (default: filter paired reads in first step)
--filterunmapped | - | Filter unmapped reads in demultiplexing step of input bam file
--krakenfilter | N | Kraken-filter with threshold N [0,1] (default: 0)
--level | N | Set BGZF compression level (default: 6)
--krakenthreads | N | Number of threads per Kraken process (default: 4)
      
### Nextflow flags
A selection of built-in Nextflow flags that may be of use (Be aware: only one `-` with these flags):
Flag | Input Type | Description
--- | --- | ---
-resume | - | Resume processing; do not re-run completed processes
-profile | NAME | Use named execution profile (see below)
-qs | N | Queue size; number of CPU cores to use (default: all)
-N | EMAIL | Send completion notification to email address
 
### Profiles   
The following execution profiles are available (use '-profile <profile>'):
Flag | Description
--- | ---
standard | execute all processes on local host with conda (default)
cluster | execute certain CPU-intensive processes on SGE cluster
docker | Run the processes inside a Docker container

# Feedback
If you have questions, feel free to write me!



