process TEST{
    label 'C4'
    container "swans:alpha"
    input:
        path(input_dir)
    output:
        stdout
    script:
    """
    ls -l $input_dir
    """

}