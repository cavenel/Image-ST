#!/usr/bin/env/ nextflow

include { QUARTONOTEBOOK as RUN_QC_DOWNSTREAM} from '../../../modules/nf-core/quartonotebook/main'

workflow QC_DOWNSTREAM {

    take:
    spatialdata               // channel: [ meta, zarr ]

    main:

    ch_versions = Channel.empty()

    //
    // Quarto reports and extension files
    //
    qc_downstream_notebook = file("${projectDir}/bin/qc_downstream.qmd", checkIfExists: true)
    extensions = Channel.fromPath("${projectDir}/assets/_extensions").collect()

    //
    // Quality controls and filtering
    //
    ch_qc_downstream_input_data = spatialdata
        .map { it -> [it[1], file("${projectDir}/bin/tissuumaps_helper.py", checkIfExists: true)] }
    ch_qc_downstream_notebook = spatialdata
        .map { it -> tuple(it[0], qc_downstream_notebook) }
    //qc_downstream_params = [
    //    input_sdata: "${meta.id}.sdata",
    //    artifact_dir: "artifacts",
    //]
    ch_qc_downstream_params = spatialdata.map { meta, zarr ->
       def params = [
            input_sdata: "${meta.id}.sdata",
            artifact_dir: "artifacts",
        ]
        return tuple(meta, params)
    }
    RUN_QC_DOWNSTREAM (
        ch_qc_downstream_notebook,
        ch_qc_downstream_params.map { it[1] },
        ch_qc_downstream_input_data,
        extensions
    )
    ch_versions = ch_versions.mix(RUN_QC_DOWNSTREAM.out.versions)
    ch_qc = RUN_QC_DOWNSTREAM.out.artifacts
        | map { meta, artifacts -> [meta, artifacts[0], meta, artifacts[1]] }
        | flatten
        | collate ( 2 )
        | branch { it ->
            tissuumaps: it[1].name.startsWith('tissuumaps')
        }
    //ch_qc_sdata = ch_qc.sdata
    ch_qc_html  = RUN_QC_DOWNSTREAM.out.html
    ch_qc_nb    = RUN_QC_DOWNSTREAM.out.notebook
    ch_qc_yml   = RUN_QC_DOWNSTREAM.out.params_yaml


    emit:
    qc_html           = ch_qc_html             // channel: [ meta, html ]
    //qc_sdata          = ch_qc_sdata          // channel: [ meta, zarr ]
    qc_nb             = ch_qc_nb               // channel: [ meta, qmd ]
    qc_params         = ch_qc_yml              // channel: [ meta, yml ]
    tissuumaps        = ch_qc.tissuumaps       // channel: [ path to folder ]

    versions          = ch_versions            // channel: [ versions.yml ]
}
