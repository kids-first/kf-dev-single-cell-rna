process UNTAR_CR {
    label 'C2'
    container 'ubuntu:latest'

    input:
    path(tar_file)

    output:
    tuple stdout, path("*")

    script:
    """
    tar xvf $tar_file \
    --ignore-failed-read \
    --wildcards "*_bc_matrix*" "*clustering*" \
    --transform  's%.*/\\([^/]*/count/.*\\)%\\1%' \
    --transform 's%count%outs%' \
    --transform 's%sample_\\([filtered|raw]\\)%\\1%g' \
    --exclude '*outs/multi*' \
    --show-transformed-names \
    | cut -f 1 -d "/" | uniq
    """
}