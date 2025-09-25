process UNTAR_REF {
    label 'C2'
    container 'ubuntu:latest'

    input:
    path(tar_file)

    output:
    path("${tar_file.name.substring(0, tar_file.name.length() - ".tar.gz".length())}")

    script:
    """
    tar xvf ${tar_file}
    """
}