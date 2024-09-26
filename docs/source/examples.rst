.. _examples-page:

Examples
========

This page provides examples for the three main ways to execute quicksand. The regular run,
a run with fixed references and a rerun with fixed references within an existing run-folder.

Please see the :ref:`quickstart-page` section to download a test-dataset (:code:`split`)
and the required datastructure (:code:`refseq`).

Regular run
~~~~~~~~~~~~

The regular run is used to get an initial overview over the taxonomic
composition of the samples. quicksand provides an overview over the detected
families and the number of ancient sequences found.

Execute quicksand like this::

    nextflow run mpieva/quicksand -r v2.2 \
        -profile   singularity \
        --split    split/ \
        --db       refseq/kraken/Mito_db_kmer22/ \
        --bedfiles refseq/genomes/ \
        --masked   refseq/masked/

The output files are grouped by family-level in the :code:`out/` directory. Extracted family-sequences
after the KrakenUniq run are stored in :code:`out/{family}/1-extracted/` while mapped, deduped and filtered sequences are saved to the
:code:`out/{family}/best/{step}/` directory after the respective processing step::

    quicksand_v2.2
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


See the :code:`final_report.tsv` for a summary of the quicksand run

Fixed references
~~~~~~~~~~~~~~~~~

quicksand is designed to work with target-enriched DNA sequences and to account for
expected families in the data. For families of interest
provide an input-file with the :code:`--fixed` flag, which specifies the reference-genomes
to use for the sequences assigned by KrakenUniq to the given family. Tags are used for the
file-names and should be unique!::

    file: fixed-references.tsv

    Taxon       Tag             Genome
    Hominidae   Homo_sapiens    /path/to/reference.fasta
    Hominidae   Another_human   /path/to/reference.fasta


and start the execution with::

    nextflow run mpieva/quicksand -r v2.2 \
        -profile   singularity \
        --split    split/ \
        --db       refseq/kraken/Mito_db_kmer22/ \
        --genomes  refseq/genomes/ \
        --bedfiles refseq/masked/
        --fixed    fixed-references.tsv

The output file structure remains the same as before. For families specified in the :code:`fixed-references.tsv` file output-files
appear in the :code:`out/{family}/fixed/{step}/` directory, together with additional output-files
that are useful in additional downstream-analyses, such as the extracted deaminated reads::

    quicksand_v2.2
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
Imagine beeing interested in the evolution of the Suidae family after having analyzed all samples with
quicksand already.

And in the final report of the analysis some lines look like this::

    Family    Species                   Reference     ReadsMapped    ProportionMapped    ReadsDeduped
    Suidae    Sus_scrofa_taivanus       best          1208           0.9028              1000

The assigned species was based on the KrakenUniq results and probably doesnt resemble the "real" species as
RefSeq contains only limited amounts of reference genomes. For any analyses that go beyond the family level, a
reanalysis with a suitable reference genome is required.

After collecting the reference genome(s) for the Suidae family, prepare a fresh fixed-references file::

    Taxon       Tag                 Genome
    Suidae      super_cool_pig      /path/to/reference.fasta

and rerun the pipeline with::

    nextflow run mpieva/quicksand -r v2.2 \
        -profile   singularity \
        --rerun    \
        --fixed    fixed-references.tsv

The (additional) output files are the ones created by the :code:`--fixed` flag::

    quicksand_v2.2
    ├── out
    │    └── {family}
    │         ├── 1-extracted
    │         │    └── {RG}_extractedReads-{family}.bam
    │         └── fixed // (family in fixed)
    │              ├── 2-aligned
    │              │     └── {RG}.{family}.{Tag}.bam
    │              ├── 3-deduped
    │              │     └── {RG}.{family}.{Tag}_deduped.bam
    │              ├── 5-deaminated
    │              │     ├── {RG}.{family}.{Tag}_deduped_deaminated_1term.bam
    │              │     └── {RG}.{family}.{Tag}_deduped_deaminated_3term.bam
    │              └── 6-mpileups
    │                    ├── {RG}.{family}.{Tag}_term1_mpiled.tsv
    │                    ├── {RG}.{family}.{Tag}_term3_mpiled.tsv
    │                    └── {RG}.{family}.{Tag}_all_mpiled.tsv
    ...
    └── final_report.tsv


The report contains now additional lines for the Suidae family with the 'fixed' references tag::

    Family    Species                   Reference     ReadsMapped    ProportionMapped    ReadsDeduped
    Suidae    Sus_scrofa_taivanus       best          1208           0.9028              1000
    Suidae    super_cool_pig            fixed         1052           0.8024              976

The final report contains a mix of best (old run) and fixed (rerun) reference entries.
