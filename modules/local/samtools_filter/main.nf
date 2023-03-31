process SAMTOOLS_FILTER {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    label 'local'
    tag "$meta.id:Flag:${params.bamfilterflag}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("filtered_*.bam"), path("filtercount.txt"), emit: bam
    path "versions.yml"                                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    samtools view -b -u $args -o filtered_${bam} ${bam}
    echo -n "\$(samtools view -c -F 128 ${bam}), \$(samtools view -c filtered_${bam})" > filtercount.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """
}