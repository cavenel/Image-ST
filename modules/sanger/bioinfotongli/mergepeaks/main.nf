process BIOINFOTONGLI_MERGEPEAKS {
    tag "$meta.id"
    label 'process_single'

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/cellgeni/tiled_spotiflow:0.5.4-1':
        'quay.io/cellgeni/tiled_spotiflow:0.5.4-1' }"

    input:
    tuple val(meta), val(ch_ind), path(csvs)

    output:
    tuple val(meta), path("${output_name}"), emit: merged_peaks
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    output_name = "${prefix}_merged_peaks_ch_${ch_ind}.wkt"
    """
    merge_overlapping_peaks.py run \
        ${csvs} \
        -output_name ${output_name} \
        ${args} \

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioinfotongli: \$(merge_overlapping_peaks.py version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ch_ind = meta.ch_ind ?: "-1"
    output_name = "${prefix}_merged_peaks_ch_${ch_ind}.wkt"
    """
    touch ${output_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioinfotongli: \$(merge_overlapping_peaks.py version)
    END_VERSIONS
    """
}
