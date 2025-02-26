#!/usr/bin/env/ nextflow

include { REGISTER_AS_SPATIALDATA } from '../subworkflows/registration'
include { MICRO_ALIGNER_REGISTRATION } from '../subworkflows/sanger/microaligner_registration/main'
include { TILED_SEGMENTATION } from '../subworkflows/sanger/tiled_segmentation/main'
include { TILED_SPOTIFLOW } from '../subworkflows/sanger/tiled_spotiflow/main' 
include { BIOINFOTONGLI_EXTRACTPEAKPROFILE as EXTRACT_PEAK_PROFILE } from '../modules/sanger/bioinfotongli/extractpeakprofile/main'  
include { POSTCODE } from '../modules/sanger/postcode/main'
include { TO_SPATIALDATA } from '../modules/local/to_spatialdata' 


workflow DECODE {
    images = channel.from(params.images)
    n_image_ch = images.map { it ->
        [it[0], it[1].size()]
    }
    MICRO_ALIGNER_REGISTRATION(images)
    TILED_SEGMENTATION(MICRO_ALIGNER_REGISTRATION.out.image, params.segmentation_method)
    TILED_SPOTIFLOW(MICRO_ALIGNER_REGISTRATION.out.image, params.chs_to_call_peaks)
    // Run the decoding
    EXTRACT_PEAK_PROFILE(MICRO_ALIGNER_REGISTRATION.out.image.join(TILED_SPOTIFLOW.out.spots_csv))
    codebook = channel.from(params.codebook).map { meta, codebook, readouts ->
        [meta,
        file(codebook, checkIfExists: true, type:'file'),
        file(readouts, checkIfExists: false, type:'file')]
    }
    POSTCODE(EXTRACT_PEAK_PROFILE.out.peak_profile.join(codebook).join(n_image_ch))
    // Contrsuct the spatial data object
    TO_SPATIALDATA(POSTCODE.out.decoded_peaks.combine(TILED_SEGMENTATION.out.wkt, by:0)
        .combine(MICRO_ALIGNER_REGISTRATION.out.image, by:0)
        .map(it ->
        [it[0], it[1], it[3], it[4]])
    )
}


workflow RNASCOPE {
    images = channel.from(params.images)
    TILED_SEGMENTATION(images, channel.from(params.segmentation_method))
    TILED_SPOTIFLOW(images, params.chs_to_call_peaks)
    TO_SPATIALDATA(TILED_SPOTIFLOW.out.spots_csv.combine(TILED_SEGMENTATION.out.wkt, by:0)
        .combine(images, by:0)
        .map(it ->
            [it[0], it[1], it[3], it[4]]
        )
    )
}