#!/usr/bin/env/ nextflow

include { REGISTER_AS_SPATIALDATA } from './subworkflows/registration'
include { MICRO_ALIGNER_REGISTRATION } from './subworkflows/sanger/microaligner_registration/main'
include { TILED_CELLPOSE } from './subworkflows/sanger/tiled_cellpose/main'
include { TILED_SPOTIFLOW } from './subworkflows/sanger/tiled_spotiflow/main' 
include { BIOINFOTONGLI_EXTRACTPEAKPROFILE as EXTRACT_PEAK_PROFILE } from './modules/sanger/bioinfotongli/extractpeakprofile/main'  
include { POSTCODE } from './modules/sanger/postcode/main'
include { TO_SPATIALDATA } from './modules/local/to_spatialdata' 

params.images = [
    [['id': "run_id"], [
        "cycle1",
        "cycle2",
    ]],
]
params.cellpose_model_dir = "/lustre/scratch126/cellgen/cellgeni/tl10/cellpose_models"
params.chs_to_call_peaks = [1, 2]
params.cell_diameters=[30]

params.codebook = [["id":'' ], "codebook.csv"]
params.chunk_size = 10000

workflow {
    images = channel.from(params.images)
    n_image_ch = images.map { it ->
        [it[0], it[1].size()]
    }
    MICRO_ALIGNER_REGISTRATION(images)
    TILED_CELLPOSE(MICRO_ALIGNER_REGISTRATION.out.image)
    TILED_SPOTIFLOW(
        MICRO_ALIGNER_REGISTRATION.out.image.combine(
            channel.from(params.chs_to_call_peaks)
        )
    )
    // Run the decoding
    EXTRACT_PEAK_PROFILE(MICRO_ALIGNER_REGISTRATION.out.image.join(TILED_SPOTIFLOW.out.spots_csv))
    if (params.codebook[2] == null) {
        readouts = null
    } else {
        readouts = file(params.codebook[2], checkIfExists: true, type:'file')
    }
    codebook = channel.from(params.codebook).map { meta, codebook, _ ->
        [meta,
        file(codebook, checkIfExists: true, type:'file'),
        readouts]
    }
    POSTCODE(EXTRACT_PEAK_PROFILE.out.peak_profile.join(codebook).join(n_image_ch))
    // Contrsuct the spatial data object
    TO_SPATIALDATA(POSTCODE.out.decoded_peaks.combine(TILED_CELLPOSE.out.wkt, by:0)
        .combine(MICRO_ALIGNER_REGISTRATION.out.image, by:0)
        .map(it ->
        [it[0], it[1], it[3], it[4]])
    )
}


workflow RNAScope {
    images = channel.from(params.images)
    TILED_CELLPOSE(images)
    TILED_SPOTIFLOW(
        images.combine(channel.from(params.chs_to_call_peaks))
    )
    TO_SPATIALDATA(TILED_SPOTIFLOW.out.spots_csv.combine(TILED_CELLPOSE.out.wkt, by:0)
        .combine(images, by:0)
        .map(it ->
            [it[0], it[1], it[3], it[4]]
        )
    )
}