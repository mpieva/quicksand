# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]
### Changed

- fix the parse_report.py script to account for cases where krakenuniq assigned a taxa without order 

## [v1.6.1] - 2023-07-04

This is a minor update to the final_report created

### Changed

- Added two columns to the end of the final_report.
  - `MeanFragmentLength`: The mean fragment length of all the DNA molecules in the bedfiltered or deduped bamfile
  - `MeanFragmentLength(3term)`: The mean fragment length of all deaminated DNA molecules in the bedfiltered or deduped bamfile
