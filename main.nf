#!/usr/bin/env/ nextflow

include { REGISTER_AS_SPATIALDATA } from './subworkflows/registration'
include { MICRO_ALIGNER_REGISTRATION } from './subworkflows/sanger/microaligner_registration/main'
include { TILED_CELLPOSE } from './subworkflows/sanger/tiled_cellpose/main'
include { TILED_SPOTIFLOW } from './subworkflows/sanger/tiled_spotiflow/main' 
include { BIOINFOTONGLI_EXTRACTPEAKPROFILE as extract_peak_profile } from './modules/sanger/bioinfotongli/extractpeakprofile/main'  
include { POSTCODE_DECODING } from './subworkflows/sanger/postcode_decoding/main' 

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
        extract_peak_profile(MICRO_ALIGNER_REGISTRATION.out.image.join(TILED_SPOTIFLOW.out.spots_csv))

        // Run the codebook conversion
        codebook = channel.from(params.codebook).map { meta, codebook, readouts ->
            [meta, file(codebook), file(readouts)]
        }
        // extract_peak_profile.out.peak_profile.join(codebook).view()
        POSTCODE_DECODING(extract_peak_profile.out.peak_profile.join(codebook))
    } else {
        transcripts = TILED_SPOTIFLOW.out.spots_csv
    }
    // transcripts.combine(TILED_CELLPOSE.out.wkt).view()
    // ASSIGN_SPOTS_TO_CELLS(TILED_SPOTIFLOW.out.spots_csv)
}