# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## WIP

### Changes

- Make the ordering in the reports consistent --> RG, Order, Family

### Bigfixes

- Fix a minor bug when re-running a quicksand analysis with the --rerun flag. The bedfilter-stats were duplicated from the respective 'best' entry and not set to '-' 
- #14: output stats for _all_ mappings done in the `stats/SAMPLE_01_mapped.tsv` and `stats/SAMPLE_02_deduped.tsv`. Only the 'best' mappings were shown before. 

## [v2.4]

### Changes

- Publish the "KrakenUniq parsed report" in the 'stats' directory
- Add an R-friendly version of the final summary report (column names w/o special characters, `R_final_report.tsv`)
- In the final_report, replace `SpeciesKmers` name with `Kmers` and include the "best" and "(family)" level stats. E.g. "4 (129)"
- for the `KmerCoverage` and `KmerDupRate` columns, also combine "best" and "(family)"
- Add a `--fixed_bedfiltering` flag to run dustmasking and bedfiltering also for fixed references (off by default)    

### Bugfixes

- Fix config-error for profile docker, causing images with ENTRYPOINT statement to fail

## [v2.3] - 2024-11-19

### Bugfixes

- Remove the default read filters in the `samtools coverage` function. (This caused results with ReadsDeduped > 0 but CoveredBP 0).
- Remove bug in the `final_report.tsv` report generation, the bug kept duplicated families (from Kraken) with MappedReads 0 in the report.
- Update default values in the report to match the format of the `final_report.tsv` (e.g. 0.0 instead of 0)  

## [v2.2] - 2024-09-26

### Changes

- Add a `--doublestranded` flag to adjust damage pattern analysis to the ones observed in data created from double stranded libraries.
  - Changes the `bam_deam_stats.py` and the `mask_deamination.py` scripts to look at 3' G to A substitutions instead of the C to T changes as done before.
  - Nothing changes for default runs

### Bugfixes 

- Fix a bug in the processing of paired fastq-files.
  - removes the "1" flag from the resulting bam files (make sure to merge your reads before quicksand!)

## [v2.1] - 2024-03-21

This version adds 4 columns to the end of the `final_report.tsv` file and adds an additional file `filtered_report_{n}p_{m}b.tsv` to the output-directory. This file should serve as a quick look on the final_report and shoult **not** be treated as the final output file

### Changes

This version alters the `final_report.tsv` file and adds an additional file `filtered_report_{n}p_{m}b.tsv` which is a filtered version of the final_report based on two freshly introduced filter-flags

- `--reportfilter_percentage` sets the filter threshold for the FamPercentage column
- `--reportfilter_breadth` sets the filter threshold for the ProportionExpectedBreadth column

Within the workflow quicksand now takes additional information from the `samtools coverage` command that analyzes the deduplicated reads. This additional information is (or is used to calculate)

- Depth of Coverage
- Breadth of Coverage
- Expected Breadth of Coverage, based on the [inStrain documentation](https://instrain.readthedocs.io/en/latest/important_concepts.html)
- Proportion of Expected Breadth, based on the [inStrain documentation](https://instrain.readthedocs.io/en/latest/important_concepts.html)

## [v2.0] - 2023-11-20

This is a rewrite of the `v1.6.1` pipeline in dsl2 syntax of nextflow
to account for nextflow-versions \>22.10

While the code was restructured, the flags, features and outputs remain the same as in `v1.6.1`
making these versions (almost) fully compatible. See the changes below.

### Changes

- instead of `FamKmers` now report `SpeciesKmers` and the respective kmer-stats to better compare assignments. Before, families with many species always had lower kmer-stats
- for each genome in the fixed-references file, run the full pipeline. They are no longer reduced to 1 reference per family as in `v1.6.1`
- parse the taxonomy directly from the DB, no need to add an additional file!
- remove the test-data and the -profile test option, as it was no longer working

## [v1.6.1] - 2023-07-04

This is a minor update to the final_report created

### Changes

- Added two columns to the end of the final_report.
  - `MeanFragmentLength`: The mean fragment length of all the DNA molecules in the bedfiltered or deduped bamfile
  - `MeanFragmentLength(3term)`: The mean fragment length of all deaminated DNA molecules in the bedfiltered or deduped bamfile
