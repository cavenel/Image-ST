VERSION="latest"

params.debug=false

process TO_SPATIALDATA {
    tag "$meta.id"
    label 'process_low'
    debug params.debug
    cache true 

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "bioinfotongli/spatialdata:${VERSION}" :
        "bioinfotongli/spatialdata:${VERSION}" }"
    publishDir params.out_dir + "/spatialdata", mode: 'copy'

    input:
    tuple val(meta), path(transcripts), path(cells_in_wkt), path(registered_image)

    output:
    tuple val(meta), path("${out_name}"), emit: spatialdata
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    out_name = "${prefix}.sdata"
    """
    export NUMBA_CACHE_DIR=/tmp/numba_cache
    /opt/conda/bin/python ${workflow.projectDir}/bin/to_spatialdata.py run \\
        --transcripts ${transcripts} \\
        --cells_in_wkt ${cells_in_wkt} \\
        --out_name ${out_name} \\
        --registered_image ${registered_image} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(to_spatialdata.py version 2>&1) | sed 's/^.*to_spatialdata.py //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sdata

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(microaligner --version 2>&1) | sed 's/^.*microaligner //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
