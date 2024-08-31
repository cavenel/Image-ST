#!/usr/bin/env/ nextflow

include { REGISTER_AS_SPATIALDATA } from './subworkflows/registration'
include { MICRO_ALIGNER_REGISTRATION } from './subworkflows/sanger/microaligner_registration/main'
include { TILED_CELLPOSE } from './subworkflows/sanger/tiled_cellpose/main'
include { TILED_SPOTIFLOW } from './subworkflows/sanger/tiled_spotiflow/main' 
include { BIOINFOTONGLI_EXTRACTPEAKPROFILE as EXTRACT_PEAK_PROFILE } from './modules/sanger/bioinfotongli/extractpeakprofile/main'  
include { POSTCODE } from './modules/sanger/postcode/main'
include { TO_SPATIALDATA } from './modules/local/to_spatialdata' 

params.debug = false
params.images = [
    [['id': "run_id"], [
        "cycle1",
        "cycle2",
    ]],
]
params.cellpose_model_dir = "/lustre/scratch126/cellgen/cellgeni/tl10/cellpose_models"
params.chs_to_call_peaks = [1, 2]
params.cell_diameters=[30]
params.decode = true
params.codebook = [["id":'' ], "codebook.csv"]
params.chunk_size = 10000


workflow {
    images = channel.from(params.images)
    MICRO_ALIGNER_REGISTRATION(images)
    TILED_CELLPOSE(MICRO_ALIGNER_REGISTRATION.out.image)
    TILED_SPOTIFLOW(
        MICRO_ALIGNER_REGISTRATION.out.image.combine(
            channel.from(params.chs_to_call_peaks)
        )
    )
    if (params.decode) {
        // Run the decoding

        //get the Peak profile
        EXTRACT_PEAK_PROFILE(MICRO_ALIGNER_REGISTRATION.out.image.join(TILED_SPOTIFLOW.out.spots_csv))

        // Run the codebook conversion
        codebook = channel.from(params.codebook).map { meta, codebook, readouts ->
            [meta,
            file(codebook, checkIfExists: true, type:'file'),
            file(readouts, type:'file')]
        }
        // extract_peak_profile.out.peak_profile.join(codebook).view()
        POSTCODE(extract_peak_profile.out.peak_profile.join(codebook))
        transcripts = POSTCODE_DECODING.out.decoded_peaks
    } else {
        transcripts = TILED_SPOTIFLOW.out.spots_csv
    }
    TO_SPATIALDATA(transcripts.combine(TILED_CELLPOSE.out.wkt, by:0)
        .combine(MICRO_ALIGNER_REGISTRATION.out.image, by:0)
        .map(it ->
        [it[0], it[1], it[3], it[4]])
    )
}