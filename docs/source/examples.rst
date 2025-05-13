.. _examples-page:

Examples
========

This page shows examples for the three main ways to execute quicksand. The regular run,
a run with fixed references and a rerun with fixed references within an existing run-folder.

Please see the :ref:`quickstart-page` section to download a test-dataset (:code:`split`)
and the required datastructure (:code:`refseq`).

Regular run
~~~~~~~~~~~~

The regular run is used to get an initial overview over the taxonomic
composition of the samples. quicksand provides an overview over the detected
families and the number of ancient sequences found.

Execute quicksand like this::

    nextflow run mpieva/quicksand -r v2.3 \
        -profile   singularity \
        --split    split/ \
        --db       refseq/kraken/Mito_db_kmer22/ \
        --bedfiles refseq/genomes/ \
        --masked   refseq/masked/

The output files are grouped by family-level in the :code:`out/` directory. Sequences binned by family-level
(after KrakenUniq) are stored in :code:`out/{family}/1-extracted/` while mapped, deduped and bedfiltered sequences are saved in the
:code:`out/{family}/best/{step}/` directories after the respective processing step::

    quicksand_v2.3
    ├── out
    │    └── {family}
    │         ├── 1-extracted
    │         │    └── {RG}_extractedReads-{family}.bam
    │         └── best
    │              ├── 2-aligned
    │              │     └── {RG}.{family}.{species}.bam
    │              ├── 3-deduped
    │              │     └── {RG}.{family}.{species}_deduped.bam
    │              └── 4-bedfiltered
    │                    └── {RG}.{family}.{species}_deduped_bedfiltered.bam
    ...
    └── final_report.tsv


See the :code:`final_report.tsv` for a summary of the quicksand run. 

Filter the final_report
~~~~~~~~~~~~~~~~~

The default quicksand-output (:code:`final_report.tsv`) is **unfiltered**, because the best 
filtering thresholds might differ between sites (and projects). However, we provide a filtered version of the report :code:`filtered_report_05p_05b.tsv`
with the default filter-thresholds applied. These thresholds are the :code:`FamPercentage` column (>=0.5%) 
and the :code:`ProportionExpectedBreadth` column (>=0.5).

Fixed references
~~~~~~~~~~~~~~~~~

quicksand is designed to work with target-enriched data. To account for
expected taxa in the sequences, users can provide a TSV-file with the :code:`--fixed` flag. This file specifies for each family the reference-genome(s) 
that quicksand uses for mapping sequences assigned by KrakenUniq to the given family. 
The 'Tags' used are used in the same way as the 'Species' (e.g. in the file-names) and should be unique!::

    file: fixed-references.tsv

    Taxon       Tag             Genome
    Hominidae   Homo_sapiens    /path/to/reference_1.fasta
    Hominidae   Another_human   /path/to/reference_2.fasta


Run quicksand with::

    nextflow run mpieva/quicksand -r v2.3 \
        -profile   singularity \
        --split    split/ \
        --db       refseq/kraken/Mito_db_kmer22/ \
        --genomes  refseq/genomes/ \
        --bedfiles refseq/masked/
        --fixed    fixed-references.tsv

The output file structure remains mostly the same. For families specified in the :code:`fixed-references.tsv` file output-files
appear in the :code:`out/{family}/fixed/{step}/` directory, together with additional output-files
that might be useful for additional downstream-analyses, such as the extracted deaminated reads::

    quicksand_v2.3
    ├── out
    │    └── {family}
    │         ├── 1-extracted
    │         │    └── {RG}_extractedReads-{family}.bam
    │         ├── best // (family not in fixed)
    |         |
    │         └── fixed // (family in fixed)
    │              ├── 2-aligned
    │              │     └── {RG}.{family}.{Tag}.bam
    │              ├── 3-deduped
    │              │     └── {RG}.{family}.{Tag}_deduped.bam
    │              ├── 4-bedfiltered #(only if --fixed_bedfiltering)
    │              │     └── {RG}.{family}.{Tag}_deduped_bedfiltered.bam
    │              ├── 5-deaminated
    │              │     ├── {RG}.{family}.{Tag}_deduped_deaminated_1term.bam
    │              │     └── {RG}.{family}.{Tag}_deduped_deaminated_3term.bam
    │              └── 6-mpileups
    │                    ├── {RG}.{family}.{Tag}_term1_mpiled.tsv
    │                    ├── {RG}.{family}.{Tag}_term3_mpiled.tsv
    │                    └── {RG}.{family}.{Tag}_all_mpiled.tsv
    ...
    └── final_report.tsv

Rerun
~~~~~~

This mode is used to repeat a run with a different set of fixed references.
For example: the final report of the analysis look like this::

    Family    Species                   Reference     ReadsMapped    ProportionMapped    ReadsDeduped
    Suidae    Sus_scrofa_taivanus       best          1208           0.9028              1000

The assigned ('best') species was based on the KrakenUniq results and might reflect the "real" species as
RefSeq contains only limited amounts of reference genomes. For any analyses that go beyond the family level, a
reanalysis with a suitable reference genome might be required.

So after collecting more reference genome(s) for the Suidae family, prepare a fresh fixed-references file::

    Taxon       Tag                 Genome
    Suidae      super_cool_pig      /path/to/reference.fasta
    Suidae      super_cool_pig2     /path/to/reference2.fasta
    Suidae      super_cool_pig3     /path/to/reference3.fasta

and rerun the pipeline with::

    nextflow run mpieva/quicksand -r v2.3 \
        -profile   singularity \
        --rerun    \
        --fixed    fixed-references.tsv

The (additional) output files are then the ones created by the :code:`--fixed` flag::

    quicksand_v2.3
    ├── out
    │    └── Suidae
    │         ├── 1-extracted
    │         │    └── {RG}_extractedReads-Suidae.bam
    │         └── fixed
    │              ├── 2-aligned
    │              │     └── {RG}.Suidae.{Tag}.bam
    │              ├── 3-deduped
    │              │     └── {RG}.Suidae.{Tag}_deduped.bam
    │              ├── 5-deaminated
    │              │     ├── {RG}.Suidae.{Tag}_deduped_deaminated_1term.bam
    │              │     └── {RG}.Suidae.{Tag}_deduped_deaminated_3term.bam
    │              └── 6-mpileups
    │                    ├── {RG}.Suidae.{Tag}_term1_mpiled.tsv
    │                    ├── {RG}.Suidae.{Tag}_term3_mpiled.tsv
    │                    └── {RG}.Suidae.{Tag}_all_mpiled.tsv
    ...
    └── final_report.tsv


The report contains now additional lines for the Suidae family with the 'fixed' references tag::

    Family    Species                   Reference     ReadsMapped    ProportionMapped    ReadsDeduped
    Suidae    Sus_scrofa_taivanus       best          1208           0.9028              1000
    Suidae    super_cool_pig            fixed         1052           0.8024              976
    Suidae    super_cool_pig2           fixed         1000           0.9001              800
    Suidae    super_cool_pig3           fixed         860            0.7551              550

The final report contains a mix of best (old run) and fixed (rerun) reference entries.
