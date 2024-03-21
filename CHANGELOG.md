# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## Unpublished

- Nothing to see here

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
