process GATHER_TAXON {
    label 'process_low'
    label 'local'
    tag "$meta.id:$meta.taxon"

    input:
    tuple val(meta), path(translate)

    output:
    tuple val(meta), path('ids.txt'), stdout, emit: ids

    script:
    """
    grep "__${meta.taxon}" ${translate} | cut -f1 | tee ids.txt | wc -l
    """
}