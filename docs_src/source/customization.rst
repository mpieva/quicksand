Customization
=============

Flags
-----

The pipeline has several required and several optional flags Required flags::

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

Optional flags::

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

Nextflow flags A selection of built-in Nextflow flags that may be of use (Be aware: only one dash - with these flags)::

    -profile  <profile>  pick a profile for executing the pipeline (see below)
    -resume              Resume processing; do not re-run completed processes
    -N        <email>    send notification to email upon fiishing
    -c        <file>     path to additional nextflow.config file for local parameters
    -w        <path>     specify the "work" directory for the nextflow intermediate files

Profiles
--------

The following profiles are available::

    To use with -profile <profile>
    standard      execute all processes with locally installed programms and conda (default, MPI-internal use only!)
    singularity   Run the pipeline in Singluarity container 
    docker        Run the pipeline in a Docker container
    cluster       execute certain CPU-intensive processes on SGE cluster


Environmental variables
-----------------------

Environmental variables can be set to reduce the number of arguments handed over to the pipeline::

    Pipeline-variables
    SED_DB        <path>     Corresponds to the --db flag
    SED_GENOME    <path>     Corresponds to the --genome flag
    SED_BEDFILES  <path>     Corresponds to the --bedfiles flag
    SED_SPECMAP   <path>     Corresponds to the --bedfiles flag

    useful nextflow-variables
    NXF_SINGULARITY_CACHEDIR <path> Where to save the pulled Singularity-images
    NXF_WORK                 <path> Corresponds to the -w flag

