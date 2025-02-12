process MERGEOUTLINES {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/cellgeni/tiled_cellpose:0.1.1':
        'quay.io/cellgeni/tiled_cellpose:0.1.1' }"

    input:
    tuple val(meta), path(outlines)

    output:
    tuple val(meta), path("${prefix}_merged.wkt"), emit: multipoly_wkts
    tuple val(meta), path("${prefix}_merged.geojson"), emit: multipoly_geojsons, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    merge_wkts.py run \\
        --sample_id ${prefix} \\
        $args \\
        $outlines

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mergeoutlines: \$(merge_wkts.py version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_merged.wkt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mergeoutlines: \$(merge_wkts.py version)
    END_VERSIONS
    """
}
