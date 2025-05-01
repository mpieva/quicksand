process RUN_DUSTMASKER{
    container (workflow.containerEngine ? "merszym/dustmasker:nextflow" : null)
    tag "${meta.Taxon}:${meta.Species}"

    input:
        tuple val(meta), path("${meta.Species}.fasta")

    output:
        tuple val(meta), path("acclist.txt"), emit: txt

    script:
        """
        dustmasker -in "${meta.Species}.fasta" -outfmt acclist > "acclist.txt"
        """
}