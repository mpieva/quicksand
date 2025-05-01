process WRITE_BEDFILES{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'quay.io/biocontainers/pandas:1.4.3' }"
    tag "${meta.Taxon}:${meta.Species}"

    input:
        tuple val(meta), path("acclist.txt")

    output:
        tuple val(meta), path("${meta.Species}.bed"), emit: bed

    script:
        """
        cat acclist.txt | \
        dustmasker_interval_to_bed.py \
        > "${meta.Species}.bed";
        """
}