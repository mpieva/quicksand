# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## Unpublished

### Bugfixes

- change in parse-report script to account for taxa without order hierarchy

## [v2.0] - 2023-11-20

This is a rewrite of the `v1.6.1` pipeline in dsl2 syntax of nextflow
to account for nextflow-versions \>22.10

While the code was restructured, the flags, features and outputs remain the same as in `v1.6.1`
making these versions fully compatible

### Changes

- instead of `FamKmers` now report `SpeciesKmers` and the respective kmer-stats to better compare assignments. Before, families with many species always had lower kmer-stats
- For each genome in the fixed-references file, run the full pipeline. **Dont** reduce to 1 reference per family as in `v1.6.1`
- parse the taxonomy directly from the DB, no need to add an additional file!
- remove the test-data

## [v1.6.1] - 2023-07-04

This is a minor update to the final_report created

### Changed

- Added two columns to the end of the final_report.
  - `MeanFragmentLength`: The mean fragment length of all the DNA molecules in the bedfiltered or deduped bamfile
  - `MeanFragmentLength(3term)`: The mean fragment length of all deaminated DNA molecules in the bedfiltered or deduped bamfile
