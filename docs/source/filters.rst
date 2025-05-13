.. _filters-page:

Filters
========

The output of quicksand is not filtered and therefore contains false-positive family-assignments. For the analysis of the quicksand-output, we recommend
applying two sets of filters.

Percentage based filter
~~~~~~~~~~~~~~~~~~~~~~~~

Filters based on a minimum percentage of sequences assigned per biological family have shown to be effective in removing misidentified sequences. 
In simulated data, we compared the relative contribution of correctly and incorrectly identified families to the total number of mapped and 
deduplicated sequences. We found that false-positive families are supported by a low percentage of sequences (median below 0.1%), 
but are generally more abundant in larger and more damaged datasets. We therefore recommend a percentage threshold of at least 0.5% of the total sequences.

The percentage of mapped and deduplicated sequences per family is listed in the :code:`FamPercentag` column of the :code:`final_report.tsv` file.


Breadth of coverage based filter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As a second filter, we calculate the evenness of coverage of mapped sequences along the reference genome for each family-level assignment. 
High evenness of coverage indicates that sequences are randomly distributed across the reference genome, which is expected when the source of the 
mapped DNA is closely related to the reference, such as from the same or a closely related species. 
In contrast, sequences from more diverged sources may only map to conserved regions, resulting in clustered alignments and low evenness of coverage. 

The parameter reported in the quicksand final summary report for evaluating coverage evenness is the 'proportion of expected breadth'. 
Breadth of coverage is defined as the proportion of the reference genome covered by at least one sequence, while Genomic coverage (or depth of coverage), 
defines the average number of times each base in the reference genome is covered by mapped sequences.

Under the assumption of random mapping to the correct reference genome, 
the breadth of coverage is a function of the genomic coverage and can be calculated using the formula empirically determined by Olm et al. 2021 [1]_.

(1) breadth of coverage = 1 - e-0.883 * coverage

We refer to the calculated breadth of coverage as the expected breadth of coverage, as it assumes mapping to the correct reference genome. 
To evaluate deviations from this expectation, we calculated for each family the proportion of expected breadth, 
defined as the ratio of the observed to the expected breadth of coverage. 
For correct family assignments, the observed breadth of coverage matches the expectations (PEB around 1), 
while the false-positive families show PEB values between 1 and 0.2. 


.. [1] Olm, M.R., Crits-Christoph, A., Bouma-Gregson, K. et al. inStrain profiles population microdiversity from metagenomic data and sensitively detects shared microbial strains. Nat Biotechnol **39**, 727–736 (2021). [https://doi.org/10.1038/s41587-020-00797-0](https://doi.org/10.1038/s41587-020-00797-0).