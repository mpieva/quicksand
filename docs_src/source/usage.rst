Usage
=====

Input
-----


Output
------

Several directories and files should have appeared upon finishing::

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
    

With the pipeline-testrun being successful, the next step is the setup of the databases

Setup databases
---------------

To run the pipeline some databases and a certain datastructure is required.

- A preindexed Kraken1-database
- mammalian mitochondrial reference genomes in the format::

    genomes/
        {family}/
            {species}.fasta(.amb/.ann ...) - fasta files preindexed with bwa
        taxid_map.tsv - A table with all nodes in the database, mapping taxid to all species within that taxon (format: '<taxid>\t<Family>\t<Species>') 
    masked:
        {species}.masked.bed - Bed files for all species in the database showing low-complexity regions

To create the structure, please run the datastructure-pipeline provided here::

    cd ..
    git clone https://github.com/MerlinSzymanski/datastructure_nf/
    nextflow run datastructure_nf/main.nf -profile singularity --outdir sediment_db 

This should create one folder within your current directory that contains all the required files ::

    sediment_db/
        kraken/
            Mito_db_kmer22/
        genomes/
            {family}/{species}.fasta
            taxid_map.tsv
        masked/
            {species}.masked.bed
            
Run the pipeline
----------------
With everything set up, the pipeline can be executed (here: directly from github)::

    nextflow run https://github.com/MerlinSzymanski/sediment_nf
         --split     <path/to/demultiplexed_reads/>
         --db        <path/to/kraken_db/>
         --genome    <path/to/reference/genomes/>
         --bedfiles  <path/to/masked_genomes/>
         --analyze
         --report
         -profile singularity

make sure that you are in a separate folder, as the output is created within the current working directory. 
see the flags section for a detailed overview of the required and optional flags!
