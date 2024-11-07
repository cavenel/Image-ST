container_version = "0.0.1"
params.debug = false

process BIOINFOTONGLI_GENERATETILES {
    tag "${meta.id}"
    debug params.debug

    label "small_mem"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "quay.io/bioinfotongli/large_image_io:${container_version}":
        "quay.io/bioinfotongli/large_image_io:${container_version}"}"

    publishDir params.out_dir + "/tile_coords"

    input:
    tuple val(meta), path(file_in)

    output:
    tuple val(meta), path("${stem}/${out_name}"), emit: tile_coords
    path "versions.yml"           , emit: versions

    script:
    stem = meta.id
    out_name = "tile_coords.csv"
    def args = task.ext.args ?: ''  
    """
    python /opt/scripts/tile_2D_image.py run \\
        --image ${file_in} \\
        --out_dir "${stem}" \\
        --out_name "${out_name}" \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(python /scripts/tile_2D_image.py version 2>&1) | sed 's/^.*tile_2D_image.py //; s/Using.*\$//' ))
    END_VERSIONS
    """
}