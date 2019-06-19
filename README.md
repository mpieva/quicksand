# Sediment pipeline (Nextflow edition)

This is a re-implentation (in Nextflow) of my original re-implementation (in
Snakemake) of the first part of Fred's sediment pipeline.

The purpose of the reimplmentation is touse Krakan for classification, and to
implement some other requested changes

## Working around two bugs in Nexflow 😔

I have stumbled across two bugs in Nextflow which affect how easy it is for you, the end user, to run a Nextflow pipeline. Neither of them affect the actual running of the pipeline itself; just the process of starting it.

I have reported both these bugs, and they are being worked on.

Until such time as they are fixed, you will have to use two small workarounds to run a Nextflow pipeline:

1. You'll have to create a conda environment that contains both Nextflow and the `conda` too, and you will have to activate this environment to run Nextflow.
2. You'll have to create a small configuration file in your home directory before you're able to run a Nextflow pipeline straight from this repository.

Let's get started:

## Preparing to run the pipeline

First, create a conda environment that contains both Nextflow and the `conda` tool itrself. Copy and paste the following command into your shell to do so. _You only have to do this once._

```bash
conda create -n nextflow -c bioconda nextflow conda
```

Now activate this new environment:

```bash
conda activate nextflow
```

You'll see that `(nextflow)` appears in front of your shell prompt when the environment is active. The environment will stay active until you close the shell session or use `conda deactivate`. The next time you log in, you will have to activate the environment again using `conda activate nextflow`.

Next, you'll need to create a configuration file `~/.nextflow/scm` that configures `vcs.eva` as a source for Nextflow pipelines. You can do so by cutting and pasting the following line into your shell. _You only have to do this once._
 
```bash
mkdir -p "${HOME}/.nextflow" && curl -s https://vcs.eva.mpg.de/visagie/sediment_nf/snippets/8/raw >"${HOME}/.nextflow/scm"
```

You can now use `nextflow pull` to pull this workflow from `vcs.eva` into your account. It is written into a special cache directory inside your home directory, so you do not have to be in any special location in the filesystem when you do this:
   
```bash
nextflow pull visagie/sediment_nf -hub eva
```

When the pipeline is updated in future, Nextflow will notify you that a newer version is available. At that point, you can update to the most recent version with the following ocmmand:
   
```bash
nextflow pull visagie/sediment_nf
```

## Running the pipeline

The pipeline requires a number of flags to be present. These are not yet fully documented.

To run the pipeline:

```bash
nextflow run visagie/sediment_nf --db </path/to/kraken.db> --bam <input bamfile> --rg <index file>
```

You can use the optional flag `--dedup` to enable deduplication of the input before Kraken is run.

If you wish to *resume* a pipeline run (e.g. if you stopped it fomr some reason, and you do not want it to redo the steps it had already completed), you need to add the flag `-resume` (note: just one `-`) after `run`:

```bash
nextflow run -resume visagie/sediment_nf --db </path/to/kraken.db> --bam <input bamfile> --rg <index file>
```
