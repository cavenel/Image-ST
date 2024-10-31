#!/usr/bin/env/ nextflow

container_version = "latest"

params.debug=false

process Spotiflow_call_peaks {
    debug params.debug
    tag "${meta.id}"

    label "gpu"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "quay.io/bioinfotongli/tiled_spotiflow:${container_version}":
        "quay.io/bioinfotongli/tiled_spotiflow:${container_version}"}"
    containerOptions = {
        workflow.containerEngine == "singularity" ? "--cleanenv --nv -B ${params.spotiflow_model_dir}:./spotiflow_models":
        ( workflow.containerEngine == "docker" ? "--gpus all -v ${params.spotiflow_model_dir}:./spotiflow_models": null )
    }
    publishDir params.out_dir + "/spotiflow_peaks"

    input:
    tuple val(meta), path(img), val(ch_ind)
    
    output:
    tuple val(meta), path("${meta.id}_ch_${ch_ind}/ch_${ch_ind}_peaks_Y*_X*.csv"), val(ch_ind), emit: peaks
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    export SPOTIFLOW_CACHE_DIR=./spotiflow_models
    /opt/conda/bin/python /scripts/Spotiflow_call_peaks.py run \
        -image_path ${img} \
        -out_dir "${meta.id}" \
        --ch_ind ${ch_ind} \
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(Spotiflow_call_peaks.py version 2>&1) | sed 's/^.*Spotiflow_call_peaks.py //; s/Using.*\$//' ))
    END_VERSIONS
    """
}


process Spotiflow_merge_peaks {
    debug params.debug
    tag "${meta.id}"

    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "quay.io/bioinfotongli/tiled_spotiflow:${container_version}":
        "quay.io/bioinfotongli/tiled_spotiflow:${container_version}"}"
    publishDir params.out_dir + "/spotiflow_peaks"

    input:
    tuple val(meta), path(csvs), val(ch_ind)

    output:
    tuple val(meta), path("${meta.id}_merged_peaks_ch_${ch_ind}.wkt"), emit: merged_peaks
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    /opt/conda/bin/python /scripts/Spotiflow_post_process.py run \
        ${csvs} \
        --ch_ind ${ch_ind} \
        --prefix "${meta.id}" \
        ${args} \

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(Spotiflow_post_process.py version 2>&1) | sed 's/^.*Spotiflow_post_process.py //; s/Using.*\$//' ))
    END_VERSIONS
    """
}


process Spotiflow_merge_channels {
    debug params.debug
    tag "${meta.id}"

    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "quay.io/bioinfotongli/tiled_spotiflow:${container_version}":
        "quay.io/bioinfotongli/tiled_spotiflow:${container_version}"}"
    publishDir params.out_dir + "/spotiflow_peaks"

    input:
    tuple val(meta), path(wkts)

    output:
    tuple val(meta), path("${meta.id}/peaks.csv"), emit: merged_channels
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    /opt/conda/bin/python /scripts/merge_wkts.py run \
        --prefix "${meta.id}" \
        ${wkts} \
        ${args} \

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(merge_wkts.py version 2>&1) | sed 's/^.*merge_wkts.py //; s/Using.*\$//' ))
    END_VERSIONS
    """
}


workflow TILED_SPOTIFLOW {
    take:
    images 

    main:
    ch_versions = Channel.empty()
    Spotiflow_call_peaks(images)
    ch_versions = ch_versions.mix(Spotiflow_call_peaks.out.versions.first())

    Spotiflow_merge_peaks(Spotiflow_call_peaks.out.peaks)
    ch_versions = ch_versions.mix(Spotiflow_merge_peaks.out.versions.first())

    Spotiflow_merge_channels(Spotiflow_merge_peaks.out.merged_peaks.groupTuple())
    ch_versions = ch_versions.mix(Spotiflow_merge_channels.out.versions.first())

    emit:
    spots_csv       = Spotiflow_merge_channels.out.merged_channels
    versions        = ch_versions
}