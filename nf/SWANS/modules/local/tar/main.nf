process TAR_OUTPUTS {
    label 'C2'
    container 'ubuntu:latest'

    input:
    tuple val(level), path(input_dir) // Use level to skip empty upstream dirs

    output:
    path("*tar.gz")

    script:
    def archive_name = level ? "${input_dir.split('/').last()}.tar.gz" : "${input_dir}.tar.gz"
    def path_name = level ?: input_dir
    """
    tar chzvf $archive_name $path_name
    """
}