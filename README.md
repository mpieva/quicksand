![MIT License](https://img.shields.io/github/license/mpieva/quicksand?style=for-the-badge)
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.11106450-ff69b4?style=for-the-badge)](https://zenodo.org/doi/10.5281/zenodo.11106450)

# quicksand

See [readthedocs](https://quicksand.readthedocs.io/en/latest/in_and_out.html) for the full documentation of the pipeline.

## Description

quicksand (**quick** analysis of **s**edimentary **an**cient **D**NA) is an open-source [Nextflow](https://doi.org/10.1038/nbt.3820) pipeline designed for rapid and accurate taxonomic classification of mammalian mitochondrial DNA (mtDNA) in aDNA samples. quicksand combines fast alignment-free classification using [KrakenUniq](https://doi.org/10.1186/s13059-018-1568-0) with downstream mapping ([BWA](https://github.com/mpieva/network-aware-bwa)), post-classification filtering, and ancient DNA authentication. quicksand is optimized for speed and portablity and requires either [Singularity](https://doi.org/10.1371/journal.pone.0177459) or [Docker](https://www.docker.com/).

## Workflow

<p align=center style="background-color:white;">
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
**Warning:** This can take several hours! For testing quicksand its recommended to just build a small database (see above)

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

Please see the [documentation](https://quicksand.readthedocs.io/en/latest/in_and_out.html) for a comprehensive description of the output files and structure!

The main summary table (`final_report.tsv`) contains one line per input file and detected family (passing the `--krakenuniq_min_kmers` and `--krakenuniq_min_reads` cutoff). The following columns are reported:

- **RG:** Name of the file analyzed
- **ReadsRaw:** Raw number of sequences in the file (paired reads only counted once)
- **ReadsFiltered** Number of sequences after filtering for `--bamfilterflag`
- **ReadsLengthfiltered:** Number of sequences after additonal filtering for sequence length (`--bamfilter_length_cutoff`)
- **Kmers:** KrakenUniq: Number of unique kmers used for classification (format: "best" and "(family)") 
- **KmerCoverage:** KrakenUniq: Kmer coverage for that classification (format: "best" and "(family)")
- **KmerDupRate:** KrakenUniq: Kmer duplication rate for that classification (format: "best" and "(family)")
- **ExtractLVL:** "f" (family) or "o" (order), set by `--taxlvl`
- **ReadsExtracted:** Number of sequences assigned by KrakenUniq
- **Order:** Detected Order
- **Family:** Detected Family
- **Species:** The reference-genome used for the mapping of 'ReadsExtracted'
- **Reference:** "best" or "fixed" (`--fixed`)
- **ReadsMapped:** Number of sequences mapped, passing the mapping quality-cutoff(`--mapbwa_quality_cutoff`)
- **ProportionMapped:** ReadsMapped / ReadsExtracted 
- **ReadsDeduped:** Number of unique sequences (removed PCR duplicates) 
- **DuplicationRate:** ReadsMapped / ReadsDeduped
- **CoveredBP:** `covbases` stat of the `samtools coverage` command. The number of bases covered in the reference genome by mapped sequences (max: ~17000)
- **ReadsBedfiltered:** Number of unique sequences after applying low-complexity bed-filtering
- **PostBedCoveredBP:** Number of bases covered in the reference genome by bedfiltered sequences
- **FamPercentage:** Percentage of unique sequences in the alignment from the total number of unique sequences in the sample (For PSF-filter)
- **Ancientness:** "-", "+" or "++". Significance level of the C-to-T or G-to-A deamination rate in the alignment (for G-to-A, use `--doublestranded`)   
- **ReadsDeam(1term):** Number of unique sequences with a C-to-T or G-to-A substitution in the _terminal_ base positions (for G-to-A, use `--doublestranded`)
- **ReadsDeam(3term):** Number of unique sequences with a C-to-T or G-to-A substitution in the _terminal three_ positions (for G-to-A, use `--doublestranded`)
- **Deam5(95ci):** The _terminal_ C-to-T substitution-rate (and the 95% confidence interval) for the 5' end
- **Deam3(95ci):** The _terminal_ C-to-T or G-to-A substitution-rate (and the 95% confidence interval) for the 3' end (for G-to-A, use `--doublestranded`)
- **Deam5Cond(95ci):** The _terminal_ C-to-T substitution-rate (and the 95% confidence interval) for the 5' end, conditioned on a substitution at the opposite end 
- **Deam3Cond(95ci):** The _terminal_ C-to-T (or G-to-A) substitution-rate (and the 95% confidence interval) for the 3' end, conditioned on a substitution at the opposite end
- **MeanFragmentLength:** Mean fragment length of all unique sequences
- **MeanFragmentLength(3term):** Mean fragment length of all 'ancient' sequences
- **Coverage:** `meandepth` of the `samtools coverage` command. Corresponds to the depth of coverage in the alignment.
- **Breadth:** `coverage` of the `samtools coverage` command / 100: The proportion of bases covered by mapped sequences in the reference genome
- **ExpectedBreadth:** Expected breadth based on 'Coverage'. Calculated using the formula: ExpectedBreadth = 1 - e^(-0.833 × Coverage)
- **ProportionExpectedBreadth:** Breadth / ExpectedBreadth (For PEB-filter)

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
> “quick analysis of sedimentary ancient DNA using quicksand”,
> [https://doi.org/10.1101/2025.08.01.668088](https://doi.org/10.1101/2025.08.01.668088).

This pipeline uses code inspired by the [nf-core](https://nf-co.re) initative, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> Ewels, Philip A., Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso, and Sven Nahnsen. 2020.
> “The Nf-core Framework for Community-curated Bioinformatics Pipelines”.
> Nature Biotechnology 38 (3): 276–78. [https://doi.org/10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x).
