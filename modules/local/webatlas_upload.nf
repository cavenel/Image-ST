process WEBATLAS_UPLOAD {
    tag "${meta.id}"
    label 'process_small'

    secret 's3_access'
    secret 's3_secret'

    input:
    tuple val(meta), path(spatialdata), val(project), val(dataset), path(config_file)

    output:
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    s3cmd cp ${config_file} s3://webatlasdataprod/${project}/${dataset}/ \
        --access_key=\${s3_access} --secret_key=\${s3_secret}

    s3cmd cp -r ${spatialdata} s3://webatlasdataprod/${project}/${dataset}/ \
        --access_key=\${s3_access} --secret_key=\${s3_secret} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        webatlas upload : \$(s3cmd --version))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sdata

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        webatlas upload : \$(s3cmd --version))
    END_VERSIONS
    """
}
