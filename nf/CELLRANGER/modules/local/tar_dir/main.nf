process TAR_DIR {
    label 'C2'
    container 'ubuntu:latest'

    input:
    tuple val(sample), path(input_dir)

    output:
    path("*tar.gz")

    script:
    """
    tar czvhf ${sample}.tar.gz ${input_dir}
    """
}