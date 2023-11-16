process GATHER_TAXON {
    label 'process_low'
    label 'local'
    tag "$meta.id:$meta.Taxon"

    input:
    tuple val(meta), path(translate)

    output:
    tuple val(meta), path('ids.txt'), stdout, emit: ids

    script:
    """
    grep "__${meta.Taxon}" ${translate} | cut -f1 | tee ids.txt | wc -l
    """
}