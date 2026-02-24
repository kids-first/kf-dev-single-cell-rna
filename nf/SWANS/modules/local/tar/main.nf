process TAR_OUTPUTS {
    label 'C2'
    container 'ubuntu:latest'

    input:
    tuple val(level), path(input_dir) // Use level to skip empty upstream dirs

    output:
    path("*tar.gz")

    script:
    // name of archive will be top level dir name or remaining dir name stripping "level" prefix if given
    archive_name = level ? "${level.split('/').last()}.tar.gz" : "${input_dir}.tar.gz"
    flags = "chzvf" // create, compress, verbose, file
    // change to path if level set
    path_name = level ? "-C ${level} ." : input_dir
    """
    tar $flags $archive_name $path_name
    """
}