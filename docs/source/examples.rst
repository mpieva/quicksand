.. _examples-page:

Examples / Tutorial
===================

On this page, we provide examples for the analysis of DNA-libraries published in three different studies to demonstrate the use of quicksand on different
library- (single- and double-stranded) and experiment types (capture and shotgun). We demonstrate how the final summary report is evaluated and how filter thresholds can be
adjusted.

The main columns used to evaluate the final summary report are `FamPercentage` and `ProportionExpectedBreadth`, which correspond to the PSF and PEB
filter values, respectively (see the `quicksand publication <https://doi.org/10.1093/molbev/msaf305>`_). 

Importantly, the PEB filter is generally more informative than the PSF filter because the breadth of coverage is a key metric for determining assignment
accuracy. However, since PEB assumes unbiased (random) mapping of sequences to the reference genome, lower values may occur when this expectation is 
violated; for instance, due to the enrichment strategy or lacking reference genome representation in the database.

Additional columns provided in the report for quality control include the KrakenUniq statistics: `Kmers` (unique kmer count), `KmerCoverage` 
(proportion of kmers observed), and `KmerDupRate` (average number of times each kmer was observed) as well as the `ProportionMapped` column. 

Download quicksand-database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please see the :ref:`quickstart-page` section to download a test-dataset (:code:`split`)
and the required datastructure (:code:`refseq`).

quicksand was developed and tested using a database built from the mammalian mtDNA genomes in RefSeq. 
For broader usability, the database provided on the FTP server was built from the full RefSeq release 221, including non-mammalian mtDNA genomes. 
This was done, because database construction is computationally intensive and may not be feasible for users with limited computational resources. 

While this expanded database enables the detection of a wider range of organisms, quicksand (and the filter thresholds) have not been 
systematically tested for the validation of non-mammalian assignments. For this section, we therefore focus on the mammalian
classifications only.

.. _ssDNA-label:

Single stranded Libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~

For the estimation of ancient DNA damage rates, quicksand assumes that sequences are derived from libraries prepared using a single stranded library 
preparation protocol (ssDNA). For the analysis of libraries prepared with a double stranded library
preparation method, see :ref:`dsDNA-label`

We first demonstrate running quicksand on ssDNA libraries enriched for mammalian mtDNA (:ref:`ssDNA-mamm`) and human mtDNA (:ref:`ssDNA-hum`) 
as well as shallow shotgun sequencing data (:ref:`ssDNA-shot`). 

As example, we use a librares created from the same sediment sample (SP3854), collected from Layer 15 of Denisova Cave, East Chamber, published in
`Slon et al. 2017 <https://doi.org/10.1126/science.aam9695>`_. The data are available in the European Nucleotide Archive (ENA), and
an overview of the libraries is provided below. An example for the analysis of multiple samples in parallel, see :ref:`dsDNA-mamm`.

+--------+-------+---------+-------------------------+-----------+------------+
| Sample | Layer | Library | Library Type            | Sequences | ENA Entry  |
+========+=======+=========+=========================+===========+============+
| SP3854 | EC 15 | R4095   | Mammalian mtDNA Capture | 1,837,960 | ERR1883545 |
| SP3854 | EC 15 | R4051   | Human mtDNA Capture     | 1,512,964 | ERR1883655 |
| SP3854 | EC 15 | R5723   | Shotgun                 | 1,638,519 | ERR1883446 |
+--------+-------+---------+-------------------------+-----------+------------+

.. _ssDNA-mamm:

ssDNA / Mammalian mtDNA Capture
*******************************

*Runtime: ~3 min.*

With the quicksand database downloaded to the :code:`refseq` directory (:ref:`quickstart-page`),
create a subdirectory for the mammalian mtDNA capture analysis and download the
required library file into a fresh input directory::

    mkdir ssDNA_mammMT && cd ssDNA_mammMT
    wget -P split
    ftp.sra.ebi.ac.uk/vol1/run/ERR188/ERR1883545/R4095.bam​

Then run quicksand with the default parameters::

    nextflow run mpieva/quicksand \
        -r v2.5 \
        --split split/ \
        --db ../refseq/kraken/Mito_db_kmer22/ \
        --genomes ../refseq/genomes/ \
        --bedfiles ../refseq/masked/ \
        -profile singularity


After completion, the final summary report is saved as :code:`quicksand_v2.5/final_report.tsv`. In total, 22 mammalian families are reported
of which 20 show :code:`Ancientness` levels of + or ++ (see Table below). Nine mammalian families pass the default PSF and PEB filter thresholds (see
:code:`quicksand_v2.5/filtered_report_0.5p_0.5b.tsv`).

+-----------------+--------------------------------+-------+-------+----------+-------+
| Family          | Best Reference                 | Reads | PSF   | Coverage | PEB   |
+=================+================================+=======+=======+==========+=======+
| Hyaenidae       | Crocuta_crocuta                | 7616  | 34.98 | 24.84x   | 0.662 |
| Rhinocerotidae  | Coelodonta_antiquitatis        | 3965  | 19.15 | 11.32x   | 0.934 |
| Bovidae         | Bos_grunniens                  | 3448  | 16.16 | 10.17x   | 0.824 |
| Ursidae         | Ursus_spelaeus                 | 1066  | 5.04  | 3.26x    | 0.757 |
| Canidae         | Canis_lupus                    | 970   | 4.52  | 2.71x    | 0.758 |
| Equidae         | Equus_ovodovi                  | 884   | 4.24  | 2.72x    | 0.703 |
| Elephantidae    | Mammuthus_columbi              | 553   | 2.65  | 1.61x    | 0.711 |
| Mustelidae      | Martes_zibellina               | 514   | 2.45  | 1.51x    | 0.795 |
| Cervidae        | Cervus_nippon                  | 281   | 1.32  | 0.82x    | 0.718 |
| Felidae         | Panthera_tigris                | 369   | 1.74  | 1.04x    | 0.346 |
| Eupleridae      | Galidia_elegans                | 225   | 1.10  | 0.7x     | 0.103 |
| Phocidae        | Monachus_monachus              | 166   | 0.81  | 0.47x    | 0.131 |
| Sciuridae       | Marmota_himalayana             | 155   | 0.75  | 0.44x    | 0.291 |
| Muridae         | Rattus_norvegicus              | 141   | 0.67  | 0.4x     | 0.193 |
| Viverridae      | Paguma_larvata                 | 112   | 0.26  | 0.32x    | 0.331 |
| Cercopithecidae | Macaca_sylvanus                | 70    | 0.33  | 0.2x     | 0.23  |
| Cricetidae      | Neodon_irene                   | 65    | 0.31  | 0.18x    | 0.443 |
| Cebidae         | Sapajus_xanthosternos          | 64    | 0.24  | 0.19x    | 0.148 |
| Hominidae       | Homo_sapiens_subsp._'Denisova' | 41    | 0.19  | 0.13x    | 1.005 |
| Otariidae       | Zalophus_japonicus             | 33    | 0.16  | 0.09x    | 0.402 |
+-----------------+--------------------------------+-------+-------+----------+-------+

The mammalian mtDNA enrichment analyzed here (targeting 242 mammalian mtDNA genomes) generally provides sufficiently random coverage across all mammalian
families for the PEB filter to be effective. The lower PEB value observed for the Hyaenidae assignment (0.66), however, is likely due to the absence of Hyaenidae
mtDNA in the capture probe set.

The Felidae, Eupleridae, Phocidae, and Sciuridae assignments pass the PSF filter but not the PEB filter. 
With the observed genomic coverages (between 0.4x and 1x), these cases are likely false positives. In contrast, the Hominidae assignment passes the PEB 
filter (PEB = 1.0) but not the PSF filter. The high PEB value for the Hominidae indicates an authentic assignment that would benefit from additional data for verification (human mtDNA capture).

In summary, for generic mammalian mtDNA capture libraries, the default quicksand filters (PSF = 0.5 and PEB = 0.5) provide a reliable set of family
assignments.

.. _ssDNA-hum:

ssDNA / Human mtDNA Capture
****************************

*Runtime: ~3 min.*

With the quicksand database downloaded to the :code:`refseq` directory (:ref:`quickstart-page`),
create a subdirectory for the mammalian mtDNA capture analysis and download the
required library file into a fresh input directory::

    mkdir ssDNA_humMT && cd ssDNA_humMT
    wget -P split ftp.sra.ebi.ac.uk/vol1/run/ERR188/ERR1883655/R4051.bam


Then run quicksand with the default parameters::

    nextflow run mpieva/quicksand \
        -r v2.5 \
        --split split/ \
        --db ../refseq/kraken/Mito_db_kmer22/ \
        --genomes ../refseq/genomes/ \
        --bedfiles ../refseq/masked/ \
        -profile singularity


After completion, the final summary report is saved as :code:`quicksand_v2.5/final_report.tsv`. A total of 15 mammalian families are
detected (Table below). Only the Hominidae family passes the default PSF and PEB filter thresholds (:code:`quicksand_v2.5/filtered_report_0.5p_0.5b.tsv`).

+-----------------+--------------------------------+-------+-------+----------+-------+
| Family          | Best Reference                 | Reads | PSF   | Coverage | PEB   |
+=================+================================+=======+=======+==========+=======+
| Hyaenidae       | Crocuta_crocuta                | 2894  | 55.97 | 9.75x    | 0.19  |
| Rhinocerotidae  | Coelodonta_antiquitatis        | 587   | 12.61 | 1.94x    | 0.192 |
| Bovidae         | Bos_taurus                     | 280   | 5.40  | 0.86x    | 0.231 |
| Equidae         | Equus_ovodovi                  | 154   | 3.24  | 0.5x     | 0.252 |
| Felidae         | Neofelis_diardi                | 146   | 2.81  | 0.43x    | 0.107 |
| Sciuridae       | Dremomys_rufigenis             | 132   | 2.60  | 0.35x    | 0.19  |
| Hominidae       | Homo_sapiens_subsp._'Denisova' | 118   | 2.28  | 0.4x     | 0.968 |
| Cervidae        | Cervus_hanglu                  | 102   | 2.17  | 0.3x     | 0.202 |
| Muridae         | Rattus_norvegicus              | 99    | 2.10  | 0.29x    | 0.171 |
| Cercopithecidae | Macaca_nigra                   | 85    | 1.80  | 0.27x    | 0.152 |
| Ursidae         | Ursus_spelaeus                 | 81    | 1.44  | 0.27x    | 0.354 |
| Mustelidae      | Martes_zibellina               | 75    | 1.50  | 0.25x    | 0.314 |
| Canidae         | Canis_aureus                   | 64    | 1.35  | 0.2x     | 0.312 |
| Elephantidae    | Mammuthus_primigenius          | 55    | 1.18  | 0.15x    | 0.322 |
| Cricetidae      | Phodopus_sungorus              | 36    | 0.77  | 0.09x    | 0.351 |
+-----------------+--------------------------------+-------+-------+----------+-------+

Even in the human mtDNA capture data, assignments are dominated by faunal “bycatch,” i.e., mtDNA from other mammalian families that is also recovered by
the human mtDNA probes. However, unlike for the mammalian mtDNA capture, non-human assignments do not meet the expectation of random mapping
because they are biased towards the subset of sequences similar to the human mtDNA probes. 

This substantially lowers the PEB-values of all non-human assignments, ranging from 0.1 to 0.35. As a consequence, the **PEB filter cannot be used
to identify true-positive assignments among the faunal bycatch**; it can do so only for the Hominidae - the only family that meets the 
criterion of random mapping.

In addition, the PSF-values in single-species captures are less informative for distinguishing true- from false-positive families. 
The above-mentioned capture bias, combined with a reduced number of total detected families, changes the distribution
of PSF values compared to the mammalian captures. We find an increased PSF for the dominant family (e.g., Hyaenidae: 55.97 vs. 34.98) and a leveling of PSF values across
the remaining families. In this example, all mammalian families pass the default PSF threshold.

Because the human mtDNA capture experiment is intentionally biased toward Hominidae, applying the PSF filter is optional, as it could remove the 
sequences we are specifically trying to detect. Instead, we recommend relying on the PEB filter to assess the quality of Hominidae assignments. 

For a more conservative approach, a stricter PEB threshold (e.g., 0.7) combined with a minimum number of reads or genomic
coverage (e.g., >0.1x) could be applied.

.. _ssDNA-shot:

ssDNA / Shotgun Sequencing
**************************

*Runtime: ~3 min.*

With the quicksand database downloaded to the :code:`refseq` directory (:ref:`quickstart-page`),
create a subdirectory for the mammalian mtDNA capture analysis and download the
required library file into a fresh input directory::

    mkdir ssDNA_shotgun && cd ssDNA_shotgun
    wget -P split ftp.sra.ebi.ac.uk/vol1/run/ERR188/ERR1883446/R5723.bam

Then run quicksand with the default parameters::

    nextflow run mpieva/quicksand \
        -r v2.5 \
        --split split/ \
        --db ../refseq/kraken/Mito_db_kmer22/ \
        --genomes ../refseq/genomes/ \
        --bedfiles ../refseq/masked/ \
        -profile singularity


After completion, the final summary report is saved as :ref:`quicksand_v2.5/final_report.tsv`. A single mammalian family, the Hyaenidea is
detected with 33 sequences (Table below), showing an Ancientness level of ++ and passing both the default PSF and PEB filter thresholds
(:code:`quicksand_v2.5/filtered_report_0.5p_0.5b.tsv`).

+-----------+-----------------+-------+-------+----------+-------+
| Family    | Best Reference  | Reads | PSF   | Coverage | PEB   |
+===========+=================+=======+=======+==========+=======+
| Hyaenidae | Crocuta_crocuta | 33    | 40.50 | 0.08x    | 1.153 |
+-----------+-----------------+-------+-------+----------+-------+

With the low sequencing depth (1,638,519 sequences), only the Hyaenidae are detected by quicksand. The assignments in the shotgun data are unaffected by
capture bias, resulting in a PEB value of ~1 for the Hyaenidae.

The final summary report includes only assignments that pass the initial KrakenUniq filters. 
As we demonstrated with simulated data (see the `quicksand publication <https://doi.org/10.1093/molbev/msaf305>`_), the default
thresholds for the number of unique kmers can be lowered for a more sensitive screening, allowing additional families to pass the initial KrakenUniq step and proceed
to mapping and downstream evaluation.

quicksand saves the KrakenUniq report for each sample, which can be examined to identify families that would be included if the filter thresholds 
were lowered. In this library, the report (:code:`quicksand_v2.5/stats/R5723.kraken.report`) lists 10 mammalian family assignments, 
of which only the Hyaenidae meets the default thresholds of 129 unique kmers and 3 sequences (Table below).

+----------+-------+--------+------------------+
| TaxReads | Kmers | Rank   | Name             |
+==========+=======+========+==================+
| 39       | 317   | family | Hyaenidae        |
| 16       | 82    | family | Bovidae          |
| 6        | 48    | family | Rhinocerotidae   |
| 15       | 30    | family | Muridae          |
| 10       | 29    | family | Vespertilionidae |
| 15       | 24    | family | Cricetidae       |
| 7        | 16    | family | Pteropodidae     |
| 6        | 12    | family | Soricidae        |
| 2        | 11    | family | Otariidae        |
| 2        | 11    | family | Ursidae          |
+----------+-------+--------+------------------+

To override the default filters and process *all* potential mammalian families, repeat the quicksand-run with the flags :code:`--krakenuniq_min_kmers 11` and
:code:`--krakenuniq_min_reads 2` to apply more permissive KrakenUniq filter thresholds::
​
    # rename the folder, to not overwrite it​
    mv quicksand_v2.5/ quicksand_v2.5.old/
    
    nextflow run mpieva/quicksand \
        -r v2.5 \
        --split split/ \
        --db ../refseq/kraken/Mito_db_kmer22/ \
        --genomes ../refseq/genomes/ \
        --bedfiles ../refseq/masked/ \
        --krakenuniq_min_kmers 11 \
        --krakenuniq_min_reads 2 \
        -profile singularity

The lower filter thresholds increase the runtime (~10 minutes) and allow 600 families to pass the KrakenUniq step. Of these, 62 assignments have at 
least one mapped sequence, and 57 pass the PEB and PSF filters. Only four *mammalian* families pass both filters (Table below). Given the low sequence
counts per family, only the Hyaenidae shows significant Ancientness levels (++).

Although the number of sequences is very low, the detected families correspond to the four most abundant ones in the sample, as determined using the mammalian capture
library.

+----------------+-------------------------+-------+-------+----------+-------+
| Family         | Best Reference          | Reads | PSF   | Coverage | PEB   |
+================+=========================+=======+=======+==========+=======+
| Hyaenidae (++) | Crocuta_crocuta         | 33    | 17.58 | 0.085x   | 1.153 |
| Rhinocerotidae | Coelodonta_antiquitatis | 5     | 2.74  | 0.01x    | 1.1   |
| Bovidae        | Procapra_przewalskii    | 2     | 1.09  | 0.0x     | 1.135 |
| Ursidae        | Ursus_spelaeus          | 1     | 0.54  | 0.0x     | 1.134 |
+----------------+-------------------------+-------+-------+----------+-------+

To summarize, while we discourage using quicksand on shallow shotgun sequencing data to draw *biological* conclusions about the taxonomic composition 
of a sediment sample, it can provide a first impression of DNA preservation and help detect high-yielding samples for follow-up capture experiments.

In contrast, Section :ref:`ds-shot` shows the application of quicksand on a deeply sequenced double-stranded library.

.. _dsDNA-label:

Double stranded Libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~

Sequences from libraries prepared with a double-stranded (dsDNA) protocol display different deamination patterns compared to those from single-stranded libraries,
showing G-to-A substitutions at the 3′ ends of DNA sequences rather than the C-to-T substitutions typical of single-stranded libraries. 

To account for these G-to-A patterns in the Ancientness assessment, quicksand needs to be run with the :code:`--doublestranded` flag

.. _dsDNA-mamm:

dsDNA / Mammalian mtDNA Capture
*******************************

*Runtime: ~1.5h*

For the analysis of double-stranded mammalian mitochondrial capture libraries, we used published data from sediment samples collected at El Mirón Cave (Cantabria,
Spain). These data were generated by `Gelabert et al. 2025 <https://doi.org/10.1038/s41467-024-55740-7>`_ and are available in the
ENA under the accession ID PRJEB74514 (Table below). The libraries were enriched for 51 mammalian mtDNA genomes (`Tejero et al. 2024 <https://doi.org/10.1016/j.heliyon.2024.e31858>`_) and sequenced to a depth of 5–10
million reads per library. In the original study, Gelabert et al. analyzed the data using *euka* to bin sequences and assess aDNA preservation, followed by a
BLAST/MEGAN-based approach to refine taxonomic classifications to the genus and species levels.

+-------------+--------+---------+-------------+--------+---------+
| ENA Entry   | Sample | Layer   | ENA Entry   | Sample | Layer   |
+=============+========+=========+=============+========+=========+
| ERR13916465 | 1      | 121     | ERR13916459 | 16     | 122/124 |
| ERR13916464 | 2      | 122/124 | ERR13916456 | 17     | 125     |
| ERR13916460 | 3      | 122/124 | ERR13916454 | 18     | 126     |
| ERR13916458 | 4      | 122/124 | ERR13916450 | 19     | 127     |
| ERR13916455 | 5      | 125     | ERR13916448 | 20     | 128     |
| ERR13916453 | 6      | 126     | ERR13916440 | 21     | 129     |
| ERR13916452 | 7      | 127     | ERR13916436 | 22     | 130     |
| ERR13916449 | 8      | 128     | ERR13916433 | 23     | 130     |
| ERR13916439 | 9      | 129     | ERR13916430 | 24     | 130     |
| ERR13916429 | 10     | 130     | ERR13916467 | 26     | 119     |
| ERR13916432 | 11     | 130     | ERR13916443 | 28     | 128     |
| ERR13916435 | 12     | 130     | ERR13916446 | 30     | 128     |
| ERR13916466 | 13     | 121     | ERR13916445 | 32     | 128     |
| ERR13916462 | 14     | 122/124 | ERR13916442 | 34     | 128     |
| ERR13916461 | 15     | 122/124 | ERR13916438 | 36     | 130     |
|             |        |         | ERR13916437 | 38     | 130     |
+-------------+--------+---------+-------------+--------+---------+

With the quicksand database downloaded to the :code:`refseq` directory (:ref:`quickstart-page`),
create a subdirectory for the mammalian mtDNA capture analysis::

    mkdir dsDNA_miron && cd dsDNA_miron


Then, download the required 32 library files into a new input directory for parallel processing::

    mkdir split && cd split \
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/038/ERR13916438/ERR13916438.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/067/ERR13916467/ERR13916467.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/037/ERR13916437/ERR13916437.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/066/ERR13916466/ERR13916466.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/082/ERR13953082/ERR13953082.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/043/ERR13916443/ERR13916443.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/046/ERR13916446/ERR13916446.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/059/ERR13916459/ERR13916459.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/035/ERR13916435/ERR13916435.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/040/ERR13916440/ERR13916440.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/052/ERR13916452/ERR13916452.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/045/ERR13916445/ERR13916445.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/030/ERR13916430/ERR13916430.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/062/ERR13916462/ERR13916462.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/053/ERR13916453/ERR13916453.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/036/ERR13916436/ERR13916436.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/032/ERR13916432/ERR13916432.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/056/ERR13916456/ERR13916456.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/061/ERR13916461/ERR13916461.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/054/ERR13916454/ERR13916454.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/049/ERR13916449/ERR13916449.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/058/ERR13916458/ERR13916458.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/055/ERR13916455/ERR13916455.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/060/ERR13916460/ERR13916460.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/029/ERR13916429/ERR13916429.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/065/ERR13916465/ERR13916465.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/064/ERR13916464/ERR13916464.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/042/ERR13916442/ERR13916442.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/033/ERR13916433/ERR13916433.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/048/ERR13916448/ERR13916448.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/050/ERR13916450/ERR13916450.bam
    wget -nc ftp.sra.ebi.ac.uk/vol1/err/ERR139/039/ERR13916439/ERR13916439.bam
    cd ../

Then, run quicksand using the default parameters and the :code:`--doublestranded` flag::

    quicksand run mpieva/quicksand \
    -r v2.5 \
    --split split \
    --genomes refseq/genomes \
    --db refseq/kraken/Mito_db_kmer22 \
    --bedfiles refseq/masked \
    --doublestranded \
    -profile singularity


After completion, the summary report contains the results for all 32 libraries. Applying the default filters (PSF of 0.5 and PEB of 0.5;
:code:`quicksand_v2.5/filtered_report_0.5p_0.5b.tsv`), recovers the same mammalian taxa reported by Gelabert et al. (Table below). 

Interestingly, we detect ancient DNA from several bird families across multiple samples: Accipitridae (11 samples), Columbidae (6), Falconidae (2), 
Phasianidae (3), and Strigidae (11), as well as ancient fish mtDNA (Salmonidae) in 14 samples, consistent with archaeological
evidence for these families in El Mirón Cave. All these families pass the PEB-filter, even though they are not part of the
mammalian enrichment set (Tejero et al. 2024; supplementary data).

+----------------+-----------------------------+------------------------------------------------+--------------------------------------------+
| Ancient Family | Best References             | Number of Samples                              | Layers                                     |
|                | (genera)                    | (Samples with > 5x Coverage)                   |                                            |
+================+=============================+================================================+============================================+
| Bovidae        | Capra, Bos                  | 31 (26)                                        | All                                        |
| Canidae        | Canis, Cuon, Vulpes         | 32 (23)                                        | All                                        |
| Cervidae       | Capreolus, Cervus, Rangifer | 32 (28)                                        | All                                        |
| Cricetidae     | Microtus                    | 28 (15)                                        | 121, 122/124, 125, 126, 127, 128, 129, 130 |
| Elephantidae   | Mammuthus                   | 4 (0)                                          | 128,130                                    |
| Equidae        | Equus                       | 22 (3)                                         | 121, 122/124, 125, 126, 127, 128, 129, 130 |
| Felidae        | Lynx, Panthera, Felis       | 24 (3)                                         | 119, 122/124, 125, 126, 127, 128, 129, 130 |
| Hominidae      | Homo                        | 13 (2)                                         | 121, 122/124, 125, 126, 127, 129           |
| Hyaenidae      | Crocuta                     | 19 (9)                                         | 119, 121, 125, 126, 128, 129, 130          |
| Leporidae      | Lepus                       | 13 (0)                                         | 121, 122/124, 125, 126, 127, 128, 129      |
| Muridae        | Apodemus                    | 5 (0)                                          | 130                                        |
| Mustelidae     | Mustela                     | 2 (2)                                          | 122/124, 128                               |
| Rhinocerotidae | Coelodonta                  | 3 (0)                                          | 130                                        |
| Soricidae      | Sorex                       | 7 (2)                                          | 119, 127, 128, 130                         |
| Suidae         | Sus                         | 3 (0)                                          | 130                                        |
| Talpidae       | Talpa                       | 13 (0)                                         | 121, 122/124, 126, 127, 128, 130           |
| Ursidae        | Ursus                       | 24 (2)                                         | 119, 122/124, 125, 126, 127, 128, 129, 130 |
+----------------+-----------------------------+------------------------------------------------+--------------------------------------------+

The original study reported 70 assignments with coverage >5x. Re-analysis of the data using quicksand identifies more samples exceeding this threshold (115) and with
higher coverage for all reported families.

To illustrate how quicksand output can be analysed further, we here repeat the human haplogroup analysis for three samples from the Solutrean levels 121, 122, and 126
reported by Gelabert et al. 2025. In all three cases, quicksand yields higher mtDNA coverage than previously reported: 4.86x, 10.71x, and 11.37x compared to 2.5x, 5.6x,
and 6.4x for levels 121, 122, and 126, respectively.

Haplogroup calling with HaploGrep 3 (`Schönherr et al., 2023 <https://doi.org/10.1093/nar/gkad284>`_) requires an mtDNA
consensus sequence based on alignment to the rCRS. However, the “best” reference genome for the Hominidae family in default
quicksand runs is often the Denisovan or Neanderthal mtDNA genome. This occurs because Neanderthals and Denisovans are classified as subspecies of modern
humans in the NCBI taxonomy (https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mod
e=tree&id=9606). Consequently, if KrakenUniq matches kmers to both the modern human and Neanderthal/Denisovan mtDNA genomes, the ‘best’ assignment is placed
at the subspecies (Neanderthal/Denisovan) level rather than at the *Homo* node.

This issue can be avoided by running quicksand with a set of predefined (“fixed”) reference genomes for selected families, either directly or 
in a second step using the :code:`--rerun` flag. 

In “rerun” mode, quicksand skips the preprocessing and KrakenUniq classification and starts from the ExtractedReads for the selected “fixed” families (e.g.,
:code:`quicksand_v2.5/out/Hominidae/1-extracted/*.bam`). quicksand then processes these sequences using the fixed reference genomes and updates the final
summary reports to include the statistics for the additional mappings. 

This example shows the combined use of the :code:`--fixed` with the :code:`--rerun` flag.

The input for the :code:`--fixed` flag is a tab-separated file specifying the family, a unique tag (or species name), and the absolute path to the mtDNA reference genome to be
used. The number of fixed references per family is not limited, allowing alignments and statistics to be generated for multiple reference genomes per family if needed.

Run quicksand with the :code:`--rerun` and :code:`--fixed` options to map all Hominidae sequences to the rCRS. The rCRS is the :code:`Homo_sapiens`
reference mtDNA genome that is part of the quicksand database

Prepare the fixed.tsv input file::

    export FAMILY=”Hominidae”​
    export TAG=”rCRS”
    export
    GENOME=”$PWD/../refseq/genomes/Hominidae/Homo_sapiens.fasta”​
    ​
    echo -e "Family\tTag\tGenome" > fixed.tsv
    echo -e $FAMILY\t$TAG\t&GENOME >> fixed.tsv

Repeat the quicksand run with the :code:`--rerun` and :code:`--fixed` options, as well as the :code:`--doublestranded` flag::

    quicksand run mpieva/quicksand \
        -r v2.5 \
        --fixed fixed.tsv \
        --rerun \
        --doublestranded \
        -profile singularity

The mapped and deduplicated bam-files can be found in the freshly created :code:`fixed` subdirectory of the family specific output-folder
(:code:`quicksand_v2.5/out/Hominidae/fixed/3-deduped/`). Construct a FASTA consensus sequence from the rCRS-mapped and deduplicated sequences 
of each sample using ANGSD v0.940 (`Korneliussen et al. 2014 <https://doi.org/10.1186/s12859-014-0356-4>`_), with the most common allele option 
(:code:`-doFasta 2` and :code:`-doCounts 1` flags), trimming the first and last base of each sequence to mitigate aDNA damage effects (:code:`-trim 1`).

Index the sample BAM-files::

    samtools index quicksand_v2.5/out/Hominidae/fixed/3-deduped/ERR13916465.Hominidae.rCRS_deduped.bam​
    samtools index quicksand_v2.5/out/Hominidae/fixed/3-deduped/ERR13916462.Hominidae.rCRS_deduped.bam
    samtools index quicksand_v2.5/out/Hominidae/fixed/3-deduped/ERR13916454.Hominidae.rCRS_deduped.bam


Create FASTA consensus sequences::

    angsd -out ERR13916465.fasta -i quicksand_v2.5/out/Hominidae/fixed/3-deduped/ERR13916465.Hominidae.rCRS_deduped.bam -doFasta 2 -doCounts 1 -trim 1
    angsd -out ERR13916462.fasta -i quicksand_v2.5/out/Hominidae/fixed/3-deduped/ERR13916462.Hominidae.rCRS_deduped.bam -doFasta 2 -doCounts 1 -trim 1
    angsd -out ERR13916454.fasta -i quicksand_v2.5/out/Hominidae/fixed/3-deduped/ERR13916454.Hominidae.rCRS_deduped.bam -doFasta 2 -doCounts 1 -trim 1


Next, use the HaploGrep 3 web tool (https://haplogrep.i-med.ac.at/) to assign human mtDNA haplogroups to each sample (Table below).

+-------------+-------+----------------------+--------------+
| Sample      | Layer | Gelabert et al. 2015 | This study   |
+=============+=======+======================+==============+
| ERR13916465 | 121   | R or U               | U            |
| ERR13916462 | 122   | U2'3'4'7'8'9         | U2'3'4'7'8'9 |
| ERR13916454 | 126   | U2'3'4'7'8'9         | U4’9         |
+-------------+-------+----------------------+--------------+

The detected haplogroups are consistent with those reported by Gelabert et al. 2025. Due to the increased amount of data recovered using quicksand, we were able to
refine the haplogroup assignments for sample ERR13916465 (U) and ERR13916454 (U4’9, a clade within U2'3'4'7'8'9).

While a complete re-evaluation of this dataset is beyond the scope of this example, it demonstrates that quicksand can be used to analyze mammalian capture
enrichment data generated from double-stranded DNA libraries using a different probeset. The results are consistent with those obtained in the original publication, and
quicksand classifies more sequences, allowing the identification of additional families and providing more data for mtDNA haplogroup assignment.

We also illustrate how quicksand output can be used for additional downstream analyses and how run options can be adjusted to meet the requirements of external
tools. While Gelabert et al. applied several additional quality-control and filtering steps to remove modern faunal contaminants (not replicated here), we show that quicksand
is suitable to provide a quick and easy starting point for subsequent in-depth, taxa-specific analyses.

.. _dsDNA-shot:

dsDNA / Shotgun Sequencing
**************************

This chapter demonstrates the analysis of a deeply sequenced dsDNA shotgun library with quicksand. 

For this example, we use data from sediment sample ‘SAT29’ from Satsurblia Cave (Georgia), described in `Gelabert et al. 2021 <https://doi.org/10.1016/j.cub.2021.06.023>`_, containing ~522 million
sequences. The BAM file available in the ENA includes raw reads processed with cutadapt v2.7 to remove sequencing adapters and poly-A–tailed reads. 

With the quicksand database downloaded to the :code:`refseq` directory (:ref:`quickstart-page`),
create a subdirectory for the shotgun library analysis and download the
required library file into a fresh input directory::

    mkdir dsDNA_sat && cd dsDNA_sat
    ​
    wget -P split \
    ftp://ftp.sra.ebi.ac.uk/vol1/err/ERR602/004/ERR6024164.bam​

Then run quicksand with the default parameters and the :code:`--doublestranded` flag::

    quicksand run mpieva/quicksand \
        -r v2.5 \
        --split split \
        --genomes ../refseq/genomes \
        --db ../refseq/kraken/Mito_db_kmer22 \
        --bedfiles ../refseq/masked \
        --doublestranded \
        -profile singularity

After completion, the final summary report is saved as :code:`quicksand_v2.5/final_report.tsv`. quicksand detects 1,237 families, of which
382 have at least one mapped sequence, and 14 pass the PSF and PEB filter thresholds. Five *mammalian* families are detected: Canidae (best reference: Canis
lupus), Bovidae (best reference: Bison bonasus), Hominidae (best reference: Homo sapiens subsp. ‘Denisova’), and Cervidae (best reference: Cervus canadensis), 
with Ancientness ratings of ++ and mtDNA coverages of 2.96x, 1.84x, 1.94x, and 0.37x, respectively. 

In addition, Cricetidae (best reference: Microtus ilaeus) is detected with Ancientness + and mtDNA coverage of 0.2x (Table below).

+------------+-----------+-------------+--------------+--------------------------------+--------+---------+-------------+-------+----------+-------+
| Family     | Extracted | Kmers       | KmerDupRate  | Best Reference                 | Mapped | Deduped | Ancientness | PSF   | Coverage | PEB   |
+============+===========+=============+==============+================================+========+=========+=============+=======+==========+=======+
| Canidae    | 2223      | 101 (8925)  | 1.69 (2.63)  | Canis_lupus_familiaris         | 1666   | 1130    | ++          | 12.42 | 2.96     | 0.950 |
| Bovidae    | 1783      | 1224 (6616) | 2.69 (2.33)  | Bison_bonasus                  | 839    | 570     | ++          | 6.18  | 1.84     | 0.989 |
| Cricetidae | 259238    | 253 (1845)  | 1.38 (613.0) | Microtus_ilaeus                | 99     | 66      | +           | 0.71  | 0.21     | 1.093 |
| Hominidae  | 1080      | 27 (7082)   | 1.3 (2.35)   | Homo_sapiens_subsp._'Denisova' | 911    | 644     | ++          | 7.10  | 1.94     | 1.001 |
| Cervidae   | 455       | 58 (1920)   | 2.21 (1.58)  | Cervus_canadensis              | 191    | 133     | ++          | 1.51  | 0.37     | 1.081 |
+------------+-----------+-------------+--------------+--------------------------------+--------+---------+-------------+-------+----------+-------+

We use this example to illustrate how the :code:`Kmers` and :code:`KmerDupRate` columns in the final summary report can be used for extended trouble-shooting. 
As shown in Table, the Cricetidae stands out with 259,238 sequences classified by KrakenUniq (Extracted). However, these sequences are based on only 
1,845 unique kmers and exhibit a very high kmer duplication rate compared to the other families (613). Particularly telling is that only 99 of the 
259,238 sequences successfully map to the best reference genome (Mapped). This combination of data indicates that sequences in the dataset containing a 
highly repetitive motif are being assigned to Cricetidae. Visual inspection of the sequences assigned to Cricetidae 
(:code:`quicksand_v2.5/out/Cricetidae/1-extracted/ERR6024164_extractedReads-Cricetidae.bam`) shows that most classified sequences are identical in length
(96 bp) and end with the same or highly similar sequence followed by a poly-G tail (ACTCCAGTCACCAGGAGGATCTCGTATGCCGTCTTCTGCTTGAAAA–PolyG), consistent with 
technical artifacts such as non-adapter-clipped reads. Although the 99 mapped Cricetidae sequences pass the PSF and PEB quicksand filters and appear 
authentic, this family-level assignment should be interpreted with caution. 

We note that Gelabert et al. (2021) used only 226,880,778 sequences for their screening, after applying additional quality-filtering steps that 
are not replicated here that would probably remove these sequences.

The assignments of the Bovidae, Canidae, and Hominidae families are consistent with the analysis by Gelabert et al. 2021, who used CENTRIFUGE in combination with the
whole non-redundant NCBI nucleotide database and a 2% classified-reads cutoff to screen the shotgun data to decide on organisms for follow-up targeted capture
experiments and whole genome mapping. In this example, we show that in this case quicksand could be used to replace the CENTRIFUGE step for the detection of mammalian
families in the sample. quicksand detects an additional family, the Cervidae, for which archaeological evidence was found in Satsurblia Cave.

In summary, this section demonstrates that quicksand can be run with default parameters on deeply sequenced shotgun data to extract the 
mammalian component that is consistent with the screening of the same data using a nuclear DNA database. 

We also show how the KrakenUniq summary statistics in the final quicksand report can assist in trouble-shooting data-anomalies.