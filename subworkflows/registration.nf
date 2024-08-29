#!/usr/bin/env/ nextflow

include { MICRO_ALIGNER_REGISTRATION } from '../subworkflows/sanger/microaligner_registration/main'

VERSION="0.5.1"
process TO_SPATIALDATA {
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "haniffalab/webatlas-pipeline:${VERSION}" :
        "haniffalab/webatlas-pipeline:${VERSION}" }"
    storeDir "${params.out_dir}/spatialdata"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path(out_name), emit: spatialdata

    script:
    out_name = "${meta.id}.sdata"
    """
    to_spatialdata.py run \
        --registered_image ${image} \
        --out_name ${out_name}
    """
}


workflow REGISTER_AS_SPATIALDATA {
    take:
    images

    main:
    MICRO_ALIGNER_REGISTRATION(images)
    TO_SPATIALDATA(MICRO_ALIGNER_REGISTRATION.out.image)

    emit:
    spatialdata = TO_SPATIALDATA.out.spatialdata
    tif         = MICRO_ALIGNER_REGISTRATION.out.image
}