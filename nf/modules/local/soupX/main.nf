process SOUPX {
    label 'C4'
    container "swans:alpha"

    input:
        val(sample)
        path(input_dir)
        val(data_type)
        val(output_dir)
        val(starting_data)
    output:
    path("${sample}_soupX/")

    script:
    """
    soupX.R \\
    $params.r_lib_path \\
    $sample \\
    $data_type \\
    $params.project \\
    $input_dir \\
    $output_dir \\
    $starting_data \\
    && mkdir ${sample}_soupX \\
    && mv $output_dir/$params.project/$sample/* ${sample}_soupX/
    """
}