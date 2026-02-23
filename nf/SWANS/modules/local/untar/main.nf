process UNTAR_CR {
    label 'C2'
    container 'ubuntu:latest'

    input:
    tuple val(meta), path(tar_file)

    output:
    tuple val(meta), path("*")

    script:
    cr_tar_args = "--ignore-failed-read " +
    "--wildcards '*_bc_matrix*' '*clustering*' " +
    "--transform  's%.*/\\([^/]*/count/.*\\)%\\1%' " +
    "--transform 's%count%outs%' " +
    "--transform 's%sample_\\([filtered|raw]\\)%\\1%g' " +
    "--exclude '*outs/multi*' " +
    "--show-transformed-names"
    tar_args = meta.input_type.contains("cellranger") ? cr_tar_args : ""
    // replace output type with dir prefix in metadata
    """
    tar xvf $tar_file \\
    $tar_args
    """
}
