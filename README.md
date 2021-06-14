# Sediment_nf
**Table of contents**
1. [Description](#description)
2. [Installation](#installation)
3. [Usage](#usage)\
3.1 [Prerequisites](#prerequisites)\
3.2 [Test](#test-the-pipeline)\
3.3 [Set up Databases](#setup-databases)\
3.4 [Running the pipeline](#run-the-pipeline)
4. [Customization](#customization)\
4.1 [Flags](#flags)\
4.2 [Profiles](#profiles)\
4.3 [Environmental variable](#environmental-variables)
5. [Pitfalls](#pitfalls)
6. [License](#license)

## Description
This repository contains an mammlian mitochondrial ancientDNA analysis pipeline.
With the overall workflow of the pipeline described in _Slon et al._, 2017, this implementation of the pipeline uses KrakenUniq for classification and modifies several steps. The pipeline is written in Nextflow and makes use of Singularity for reproducibility.

<p align=center>
    <img src="assets/pipeline_overview.png" alt="Graphical overview over the pipeline workflow" width=500 height=500>
</p>
    
The **Input** of the pipeline consists of demultiplexed, trimmed and overlap-merged reads in either .bam or .fastq format. For each readgroup, the pipeline **outputs** the raw reads in .bam-format at each stage of the pipeline. Additionally the number of reads extracted, mapped, deduplicated, and bedfiltered as well as the estimated DNA-damage is reported for a quick overview in one big summary-file. See the [Test](#test-the-pipeline) section for an overview over the output-files. Taxonomic assignments are evaluated on a family level 

## Installation

To download the pipeline, please open your Terminal and type:
```git clone github.com/MerlinSzymanski/sediment_nf```

## Usage
### Prerequisites

To run the pipeline you need the two following programms installed
1. Nextflow (tested on v20.04.10): [Installation](https://www.nextflow.io/docs/latest/getstarted.html)
2. Singularity (tested on v3.7.1): [Installation](https://sylabs.io/guides/3.0/user-guide/installation.html)\
or Docker: [Installation](https://docs.docker.com/get-docker/)

Singularity/Docker use containers to run software in, which ensures an consistent environment. The image that gets pulled and used by the pipeline is hosted on Dockerhub and can be inspected [here](https://hub.docker.com/r/merszym/sediment_nf)

### Test the pipeline
If  `nextflow` and `docker or singularity` are installed, on can **test** the pipeline. To do that, unpack the test database provided in the repository
```bash
mkdir test_run && cd test_run
tar -xvzf ../sediment_nf/assets/test/kraken/database.tar.gz
```
and run the test set
```
nextflow run ../sediment_nf/main.nf \
    --split ../sediment_nf/assets/test/split/ \
    --genome ../sediment_nf/assets/test/genomes/ \
    --bedfiles ../sediment_nf/assets/test/masked/ \
    --db TestDB/ \
    --specmap ../sediment_nf/assets/test/genomes/specmap.tsv \
    --analyze \
    --report \
    -profile singularity
```
exchange `-profile singularity` with `-profile docker` if the choice fell on docker over singularity. Also be aware that the `-profile` flag has only one dash!
To have an overview over the flags, see [Flags](#flags)

Several directories and files should have appeared upon finishing:
```
out/
    {family}/
        {readgroup}_extractedReads-{Family}.bam - All reads assigned by KrakenUniq to that family
        aligned/
            {readgroup}.{species}.bam - All extractedReads from that family mapped against the genome of that species
            {readgroup}.{species}_deduped.bam - the {readgroup}.{species}.bam file, but depleted of PCR duplicates
        bed/
            {readgroup}.{species}_deduped_bedfiltered.bam - the {readgroup}.{species}_deduped.bam file but additionally depleted of reads overlapping low-complexity regions
kraken/
    {readgroup}.report - the raw krakenUniq report
    {readgroup}.translate - the translated krakenUniq report (MPA-Format)
stats/
    splitcounts.tsv - contains the number of reads per readgroup
    {readgroup}_extracted.tsv - contains the number of reads per extracted family
    {readgroup}_mapped.tsv - contains the number of reads mapped against the reference genome
    {readgroup}_mapped_coverage.tsv - contains the number of covered basepairs of each mapping
    {readgroup}_unique_mapped.tsv - contains the number of deduplicated reads mapped against the reference genome
    {readgroup}_bedfiltered.tsv - contains the number of reads that remain in the bam-file after bedfiltering of low-complexity regions
    {readgroup}_deamination_stats - contain an estimation of the 'ancientness' of families based on deamination frequencies of recovered reads
reports/ - contains stats about the nextflow run
    report.html
    timeline.html
    trace.tsv
work/ - contains intermediate files required by nextflow. Can be deleted after the run has finished
final_report.tsv - a summary of all the stats in stats/ for families that made it past the mapping as an easy to parse tsv-file
```
With the pipeline-testrun being successful, the next step is the **setup of the databases**

### Setup databases

To run the pipeline some databases and a certain datastructure is required.

```
- A preindexed Kraken1-database
- mammalian mitochondrial reference genomes in the format:
    genomes/
        {family}/
            {species}.fasta(.amb/.ann ...) - fasta files preindexed with bwa
        taxid_map.tsv - A table with all nodes in the database, mapping taxid to all species within that taxon (format: '<taxid>\t<Family>\t<Species>') 
    masked:
        {species}.masked.bed - Bed files for all species in the database showing low-complexity regions
```
To create the structure, please run the datastructure-pipeline provided [here](https://github.com/MerlinSzymanski/datastructure_nf/)
```
cd ..
git clone https://github.com/MerlinSzymanski/datastructure_nf/
nextflow run datastructure_nf/main.nf -profile singularity --outdir sediment_db 
```
This should create one folder within your current directory that contains all the required files
```
sediment_db/
    kraken/
        Mito_db_kmer22/
    genomes/
        {family}/{species}.fasta
        taxid_map.tsv
    masked/
        {species}.masked.bed
```

### Run the pipeline

With everything set up, the pipeline can be executed (here: directly from github):
```
nextflow run https://github.com/MerlinSzymanski/sediment_nf
     --split     <path/to/demultiplexed_reads/>
     --db        <path/to/kraken_db/>
     --genome    <path/to/reference/genomes/>
     --bedfiles  <path/to/masked_genomes/>
     --analyze
     --report
     -profile singularity
```
make sure that you are in a separate folder, as the output is created within the current working directory.
see the [flags](#flags) section for a detailed overview of the required and optional flags!

## Customization

### Flags:
The pipeline has several required and several optional flags
**Required flags**
```
     --db        <path/to/kraken_db/>
     --genome    <path/to/reference/genomes/>
     --bedfiles  <path/to/masked_genomes/>
     
     EITHER
     --split     <path/to/demultiplexed_reads/> 
                 Directory with already demultiplexed, overlap-merged and trimmend bam- or fastqfiles. 
                 Format 'split/{readgroup}.fastq'
     OR
     --bam       <path/to/still/multiplexed/bamfile.bam>
     --rg        <path/to/list/of/used/indices.tsv>
                 Format 'Readgroup\tP7\tP5'
```
**Optional flags**
```
--analyze             include deamination stats and reports
--report              create the final_report.tsv
--specmap        <path/to/specmap.tsv>
                      Ignore the "best species" assigned by the pipeline. Always map extractedReads to the species specified here
                      Format: Family\tSpecies_one,Species_two
                      Species_name must correspond to the filename in the genomes diretory
--capture        <Family,Family>
                      Apply reduced filters to the families (skip Bedfiltering)
--cutoff         <N>  Length cutoff after BWA-mapping (default: 35)
--quality        <N>  Mapping quality filter after BWA-mapping (default:25)
--min_kmers      <N>  Minimum required unique kmers assigned by KrakenUniq to a family to keep it (default: 129)
--min_reads      <N>  Minimum required assigned reads by krakenUniq to a family to keep it (default: 3)
--keeppaired          Dont filter out unmerged reads (default: filterpaired)
--filterunmapped      Filter out unmapped reads (during demultiplexing premapped --bam input file)
--level          <N>  Set BGZF compression level (default: 6)
--krakenthreads  <N>  Number of threads per Kraken process (default: 4)
```
**Nextflow flags**
A selection of built-in Nextflow flags that may be of use (Be aware: only one dash `-` with these flags):
```
-profile  <profile>  pick a profile for executing the pipeline (see below)
-resume              Resume processing; do not re-run completed processes
-N        <email>    send notification to email upon fiishing
-c        <file>     path to additional nextflow.config file for local parameters
-w        <path>     specify the "work" directory for the nextflow intermediate files
```

### Profiles
The following profiles are available 
```
To use with -profile <profile>
standard      execute all processes with locally installed programms and conda (default, MPI-internal use only!)
singularity   Run the pipeline in Singluarity container 
docker        Run the pipeline in a Docker container
cluster       execute certain CPU-intensive processes on SGE cluster
```

### Environmental variables
Environmental variables can be set to reduce the number of arguments handed over to the pipeline

```
Pipeline-variables
SED_DB        <path>     Corresponds to the --db flag
SED_GENOME    <path>     Corresponds to the --genome flag
SED_BEDFILES  <path>     Corresponds to the --bedfiles flag
SED_SPECMAP   <path>     Corresponds to the --bedfiles flag

useful nextflow-variables
NXF_SINGULARITY_CACHEDIR <path> Where to save the pulled Singularity-images
NXF_WORK                 <path> Corresponds to the -w flag
```
## Pitfalls
**Path does not exist in Singularity**: Using a shared file-system, discs might need to get bound into the container before accessing files stored on that disc. One can do that by adding a singularity-command to a custom config-file:
```
singularity {
  runOptions = "--bind /mnt"
}
``` 
and including that file in your run
```
nextflow run ... -c custom_config_file
```

## License
