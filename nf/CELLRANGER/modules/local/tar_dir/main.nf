process TAR_DIR {
    label 'C2'
    container 'ubuntu:latest'

    input:
    tuple val(prefix), path(input_dir)

    output:
    path("*tar.gz")

    script:
    """
    tar czvhf ${prefix}.tar.gz ${input_dir}
    """
}