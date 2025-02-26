process BIOINFOTONGLI_GENERATETILES {
    tag "${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "quay.io/bioinfotongli/large_image_io:0.0.2":
        "quay.io/bioinfotongli/large_image_io:0.0.2"}"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("${output_name}"), emit: tile_coords
    path "versions.yml"           , emit: versions

    script:
    stem = meta.id
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    output_name = "${prefix}_tile_coords.csv"
    """
    tile_2D_image.py run \\
        --image ${image} \\
        --output_name "${output_name}" \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioinfotongli: \$(tile_2D_image.py version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    output_name = "${prefix}_tile_coords.csv"
    """
    touch "${output_name}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioinfotongli: \$(tile_2D_image.py version)
    END_VERSIONS
    """
}