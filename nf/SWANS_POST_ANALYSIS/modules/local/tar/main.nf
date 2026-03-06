process TAR_OUTPUTS {
    label 'C2'
    container 'ubuntu:latest'

    input:
    tuple val(level), val(archive_name), path(input_dir) // Use level to skip empty upstream dirs

    output:
    path("*tar.gz")

    script:
    // name of archive will be explicitly given, be the top level dir name or remaining dir name stripping "level" prefix if given
    def tar_name = archive_name ?: level ? "${level.split('/').last()}.tar.gz" : "${input_dir}.tar.gz"
    flags = "chzvf" // create, compress, verbose, file
    // change to path if level set
    path_name = level ? "-C ${level} ." : input_dir
    """
    tar $flags $tar_name $path_name
    """
}