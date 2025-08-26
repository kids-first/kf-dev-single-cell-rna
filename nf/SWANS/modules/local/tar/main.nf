process UNTAR_CR {
    label 'C2'
    container 'ubuntu:latest'

    input:
    tuple val(sample_id), path(tar_file)

    output:
    tuple val(sample_id), path("cell_ranger_${sample_id}")

    script:
    """
    mkdir cell_ranger_${sample_id}

    tar xvf ${tar_file} -C cell_ranger_${sample_id}
    """
}   