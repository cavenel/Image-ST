container_version = "0.5.1"

params.debug = false

process BIOINFOTONGLI_TILEDSPOTIFLOW {
    debug params.debug
    tag "${meta.id}"

    label "gpu"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "quay.io/bioinfotongli/tiled_spotiflow:${container_version}":
        "quay.io/bioinfotongli/tiled_spotiflow:${container_version}"}"
    containerOptions = {
        workflow.containerEngine == "singularity" ? "--cleanenv --nv":
        ( workflow.containerEngine == "docker" ? "--gpus all": null )
    }
    publishDir params.out_dir + "/spotiflow_peaks"

    input:
    tuple val(meta), val(x_min), val(y_min), val(x_max), val(y_max), path(image), val(ch_ind)
    
    output:
    tuple val(meta), val(ch_ind), path("${out_dir}/${out_name}"), emit: peaks
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    out_dir = "${meta.id}_ch_${ch_ind}"
    out_name = "ch_${ch_ind}_peaks_Y${y_min}_X${x_min}.csv"
    """
    /opt/conda/bin/python /scripts/Spotiflow_call_peaks.py run \
        -image_path ${image} \
        -out_dir "${out_dir}" \
        -x_min ${x_min} \
        -y_min ${y_min} \
        -x_max ${x_max} \
        -y_max ${y_max} \
        -ch_ind ${ch_ind} \
        -out_name "${out_name}" \
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(/opt/conda/bin/python /scripts/Spotiflow_call_peaks.py version 2>&1) | sed 's/^.*Spotiflow_call_peaks.py //; s/Using.*\$//' ))
    END_VERSIONS
    """
}