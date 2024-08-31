params.debug=false
params.chunk_size=10000


process Codebook_conversion {
    debug params.debug
    tag "${meta.id}"

    cpus 1
    memory 100.MB

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.gmm_sif:
        'bioinfotongli/decoding:latest'}"

    storeDir params.out_dir + "/codebook_metadata"

    input:
    tuple file(codebook), val(channel_map), val(sep)

    output:
    path "taglist.csv", emit: taglist_name
    path "channel_info.csv", emit: channel_info_name
    /*path "channel_info.pickle", emit: channel_infos*/
    path "versions.yml"           , emit: versions

    script:
    """
    /scripts/codebook_convert.py -csv_file ${codebook} -channel_map "${channel_map}" -sep "${sep}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(/scripts/codebook_convert.py version 2>&1) | sed 's/^.*/scripts/codebook_convert.py //; s/Using.*\$//' ))
    END_VERSIONS
    """
}


process Get_meatdata {
    debug params.debug
    tag "${meta.id}"

    cpus 1
    memory 100.MB

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        params.gmm_sif:
        'bioinfotongli/decoding:latest'}"

    storeDir params.out_dir + "/codebook_metadata"
    /*publishDir params.out_dir + "/decoding_metadata", mode:"copy"*/

    input:
    path(taglist_name)
    path(channel_info_name)

    output:
    path "barcodes_01.npy", emit: barcodes
    path "gene_names.npy", emit: gene_names
    path "channel_info.pickle", emit: channel_infos

    path "versions.yml"           , emit: versions

    script:
    """
    /scripts/get_metadata.py -auxillary_file_dir ./  -taglist_name ${taglist_name} -channel_info_name ${channel_info_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(/scripts/get_metadata.py version 2>&1) | sed 's/^.*/scripts/get_metadata.py //; s/Using.*\$//' ))
    END_VERSIONS
    """
}


process Decode_peaks_old {
    debug params.debug
    tag "${meta.id}"
    cache true

    cpus 1

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'bioinfotongli/postcode:latest':
        'bioinfotongli/postcode:latest'}"

    publishDir params.out_dir + "/decoded", mode:"copy"

    input:
    tuple val(stem), file(spot_profile), file(spot_loc), file(barcodes_f), file(gene_names_f), file(channel_info_f)
    val(chunk_size)

    output:
    tuple val(stem), path("${stem}_decoded_df.tsv"), emit:decoded_peaks
    path "${stem}_decode_out_parameters.pickle" optional true

    path "versions.yml"           , emit: versions

    script:
    """
    /scripts/decode.py --spot_profile ${spot_profile} --spot_loc ${spot_loc} --barcodes_01 ${barcodes_f} --gene_names ${gene_names_f} --channels_info ${channel_info_f} --stem ${stem} --chunk_size ${chunk_size}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(/scripts/decode.py version 2>&1) | sed 's/^.*/scripts/decode.py //; s/Using.*\$//' ))
    END_VERSIONS
    """
}


process Decode_peaks {
    debug params.debug
    tag "${meta.id}"
    cache true

    cpus 1

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'bioinfotongli/postcode:latest':
        'bioinfotongli/postcode:latest'}"
    publishDir params.out_dir + "/decoded"

    input:
    tuple val(meta), file(spot_profile), file(spot_loc), file(codebook), file(readouts)
    val(chunk_size)

    output:
    tuple val(meta), path("${out_name}"), emit:decoded_peaks
    tuple val(meta), path("${prefix}_decode_out_parameters.pickle"), optional: true

    path "versions.yml"           , emit: versions

    script:
    prefix = meta.id ?: "none_decoded"
    out_name = "${prefix}_decoded_spots.csv"
    def args = task.ext.args ?: ""
    """
    /scripts/decode.py run --spot_profile_p ${spot_profile} --spot_locations_p ${spot_loc} \\
        --codebook_p ${codebook} --out_name ${out_name} --readouts_csv ${readouts} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(/scripts/decode.py version 2>&1) | sed 's/^.*/scripts/decode.py //; s/Using.*\$//' ))
    END_VERSIONS
    """
}


workflow POSTCODE_DECODING_OLD {

    take:
    codebook // channel: [ file[ codebook ], val(channel_map), val(codebook_sep) ]
    peak_profile // channel: [ val(meta), file[ peak_profile ] ]

    main:
    ch_versions = Channel.empty()
    Codebook_conversion(codebook)
    ch_versions = ch_versions.mix(Codebook_conversion.out.versions.first())

    Get_meatdata(Codebook_conversion.out.taglist_name, Codebook_conversion.out.channel_info_name)
    ch_versions = ch_versions.mix(Get_meatdata.out.versions.first())

    for_decoding = peak_profile 
        .combine(Get_meatdata.out.barcodes)
        .combine(Get_meatdata.out.gene_names)
        .combine(Get_meatdata.out.channel_infos)

    Decode_peaks(for_decoding, params.chunk_size)
    ch_versions = ch_versions.mix(Decode_peaks.out.versions.first())

    emit:
    decoded_peaks      = Decode_peaks.out.decoded_peaks          // channel: [ val(meta), file[ decoded_peaks ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

workflow POSTCODE_DECODING {

    take:
    to_decode  // channel: [ val(meta), file[ peak_profile ], file[ spot_loc], file[ codebook ] ]

    main:
    ch_versions = Channel.empty()

    Decode_peaks(to_decode, params.chunk_size)
    ch_versions = ch_versions.mix(Decode_peaks.out.versions.first())

    emit:
    decoded_peaks      = Decode_peaks.out.decoded_peaks          // channel: [ val(meta), file[ decoded_peaks ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}