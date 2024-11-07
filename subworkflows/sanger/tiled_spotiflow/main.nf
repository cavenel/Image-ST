#!/usr/bin/env/ nextflow
include { BIOINFOTONGLI_GENERATETILES as GENERATE_TILE_COORDS } from '../../../modules/sanger/bioinfotongli/generatetiles/main'
include { BIOINFOTONGLI_TILEDSPOTIFLOW } from '../../../modules/sanger/bioinfotongli/tiledspotiflow/main'

container_version = "0.1.0"

params.debug=false
params.chs_to_call_peaks = [2]


process Spotiflow_merge_tiled_peaks {
    debug params.debug
    tag "${meta.id}"

    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "quay.io/bioinfotongli/tiled_spotiflow:${container_version}":
        "quay.io/bioinfotongli/tiled_spotiflow:${container_version}"}"
    publishDir params.out_dir + "/spotiflow_peaks"

    input:
    tuple val(meta), val(ch_ind), path(csvs)

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
    chs_to_call_peaks

    main:
    ch_versions = Channel.empty()
    GENERATE_TILE_COORDS(images)
    images_tiles = GENERATE_TILE_COORDS.out.tile_coords.splitCsv(header:true, sep:",").map{ meta, coords ->
        [meta, coords.X_MIN, coords.Y_MIN, coords.X_MAX, coords.Y_MAX]
    }

    BIOINFOTONGLI_TILEDSPOTIFLOW(images_tiles.combine(images, by:0).combine(chs_to_call_peaks))
    ch_versions = ch_versions.mix(BIOINFOTONGLI_TILEDSPOTIFLOW.out.versions.first())

    Spotiflow_merge_tiled_peaks(BIOINFOTONGLI_TILEDSPOTIFLOW.out.peaks.groupTuple(by:[0,1]))
    ch_versions = ch_versions.mix(Spotiflow_merge_tiled_peaks.out.versions.first())

    Spotiflow_merge_channels(Spotiflow_merge_tiled_peaks.out.merged_peaks.groupTuple())
    ch_versions = ch_versions.mix(Spotiflow_merge_channels.out.versions.first())

    emit:
    spots_csv       = Spotiflow_merge_channels.out.merged_channels
    versions        = ch_versions
}