include { BIOINFOTONGLI_CELLPOSE as CELLPOSE } from '../../../modules/sanger/bioinfotongli/cellpose/main'
include { BIOINFOTONGLI_GENERATETILES as GENERATE_TILE_COORDS } from '../../../modules/sanger/bioinfotongli/generatetiles/main'

params.images = [
    [["id":"test1"], "file1"],
    [["id":"test2"], "file2"],
]
params.cell_diameters = [30, 40]
params.cellpose_model_dir = "/lustre/scratch126/cellgen/cellgeni/tl10/cellpose_models"
params.debug=false

container_version = "0.1.0"


process MERGE_OUTLINES {
    tag "${meta.id}"
    debug params.debug

    label "medium_mem"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "quay.io/bioinfotongli/tiled_cellpose:${container_version}":
        "quay.io/bioinfotongli/tiled_cellpose:${container_version}"}"

    publishDir params.out_dir + "/cellpose_segmentation_merged_wkt"

    input:
    tuple val(meta), val(cell_diameter), path(wkts)

    output:
    tuple val(meta), val(cell_diameter), path("${stem}_merged.wkt"), emit: merged_wkt
    tuple val(meta), val(cell_diameter), path("${stem}_merged.geojson"), emit: merged_geojson
    path "versions.yml"           , emit: versions

    script:
    stem = "${meta.id}_diam-${cell_diameter}"
    def args = task.ext.args ?: ''  
    """
    /opt/conda/bin/python /scripts/merge_wkts.py run \
        --sample_id "${stem}" \
        ${wkts} \
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(/scripts/merge_wkts.py version 2>&1) | sed 's/^.*merge_wkts.py //; s/Using.*\$//' ))
    END_VERSIONS
    """
}


workflow TILED_SEGMENTATION {
    take:
    images 
    method

    main:
    ch_versions = Channel.empty()
    GENERATE_TILE_COORDS(images)
    ch_versions = ch_versions.mix(GENERATE_TILE_COORDS.out.versions.first())

    images_tiles = GENERATE_TILE_COORDS.out.tile_coords.splitCsv(header:true, sep:",").map{ meta, coords ->
        [meta, coords.X_MIN, coords.Y_MIN, coords.X_MAX, coords.Y_MAX]
    }

    if (method == "CELLPOSE") {
        CELLPOSE(images_tiles.combine(images, by:0).combine(channel.from(params.cell_diameters)))
        wkts = CELLPOSE.out.wkts.groupTuple(by:[0,1])
        ch_versions = ch_versions.mix(CELLPOSE.out.versions.first())
    } else if (method == "STARDIST") {
        STARDIST(images_tiles)
        wkts = STARDIST.out.wkts.groupTuple(by:[0,1])
        ch_versions = ch_versions.mix(STARDIST.out.versions.first())
    } else {
        error "Invalid segmentation method: ${method}"
    }
    MERGE_OUTLINES(wkts)
    ch_versions = ch_versions.mix(MERGE_OUTLINES.out.versions.first())

    emit:
    wkt         = MERGE_OUTLINES.out.merged_wkt   // channel: [ val(meta), [ wkt ] ]
    versions    = ch_versions                     // channel: [ versions.yml ]
}