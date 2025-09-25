process UNTAR_CR {
    label 'C2'
    container 'ubuntu:latest'

    input:
    tuple val(sample_id), path(tar_file)

    output:
    tuple val(sample_id), path("*")

    script:
    """
    transform_cr_tar.sh ${tar_file}
    """
}