process UNTAR_CR {
    label 'C2'
    container 'ubuntu:latest'

    input:
    tuple val(source), path(tar_file)

    output:
    tuple val(source), stdout, path("*")

    script:
    def cr_tar_args = "--ignore-failed-read " +
    "--wildcards '*_bc_matrix*' '*clustering*' " +
    "--transform  's%.*/\\([^/]*/count/.*\\)%\\1%' " +
    "--transform 's%count%outs%' " +
    "--transform 's%sample_\\([filtered|raw]\\)%\\1%g' " +
    "--exclude '*outs/multi*' " +
    "--show-transformed-names"
    def tar_args = source == "cellranger" ? cr_tar_args : ""
    """
    tar xvf $tar_file \\
    $tar_args \\
    | cut -f 1 -d "/" | uniq
    """
}