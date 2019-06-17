# Sediment pipeline (Nextflow edition)

This is a re-implentation (in Nextflow) of my original re-implementation (in
Snakemake) of the first part of Fred's sediment pipeline.

The purpose of the reimplmentation is touse Krakan for classification, and to
implement some other requested changes

## Running the pipeline

### Work around a bug in Nexflow 😔

In theory you should be able to run this pipeline by referring to its URL:
```bash
nextflow run https://vcs.eva.mpg.de/visagie/sediment_nf
```
…**but this does not work** due to a bug in Nextflow. I have reported the bug, but until it is fixed the following workaround can be used:

### Installation procedure
 
1. Create a tiny file `~/.nextflow/scm` that will define a configuration profile
   for `vcs.eva`. You can do this by copying the following line and pasting it
   into your shell:
   ```bash
   mkdir -p "${HOME}/.nextflow" && curl -s https://vcs.eva.mpg.de/visagie/sediment_nf/snippets/8/raw >"${HOME}/.nextflow/scm"
   ```
   >You only have to do this once!
   
2. Pull the pipeline to your local account with the following command:
   ```bash
   nextflow pull visagie/sediment_nf -hub eva
   ```
   >You only have to do this once as well!
   
3. Refer to the pipeline as `visagie/sediment_nf` in all Nextflow commands, e.g.:
   ```bash
   nextflow run visagie/sediment_nf --db /some/database
   ```
 
4. When Nextflow notifies you that there is a newer version of the pipeline, you may update your local copy with:

   ```bash
   nextflow pull visagie/sediment_nf
   ```
   
