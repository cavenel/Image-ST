VERSION="v1.0.6"

params.debug=false

process BIOINFOTONGLI_MICROALIGNER {
    tag "$meta.id"
    label 'process_large'
    debug params.debug
    cache true 

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "bioinfotongli/microaligner:${VERSION}" :
        "bioinfotongli/microaligner:${VERSION}" }"
    publishDir params.out_dir + "/registered_stacks", mode: 'copy'

    input:
    tuple val(meta), path(config), path(images)

    output:
    tuple val(meta), path("${prefix}*_reg_result_stack.tif"), emit: registered_image
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    microaligner --config ${config} \
        --NumberOfWorkers ${task.cpus} \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(microaligner --version 2>&1) | sed 's/^.*microaligner //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_optflow_reg_result_stack.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(microaligner --version 2>&1) | sed 's/^.*microaligner //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
