process TAR_OUTPUTS {
    label 'C2'
    container 'ubuntu:latest'

    input:
    tuple val(level), path(input_dir) // Use level to skip empty upstream dirs

    output:
    path("*tar.gz")

    script:
    def archive_name = level ? "${level.split('/').last()}.tar.gz" : "${input_dir}.tar.gz"
    def flags = "chzvf" // create, compress, verbose, file
    def path_name = level ? "-C ${level} ." : input_dir
    """
    tar $flags $archive_name $path_name
    """
}