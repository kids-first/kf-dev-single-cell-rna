process SOUPX {
    label 'C4'
    container "swans:alpha"

    input:
        tuple val(sample), path(input_dir)
        val(data_type)
        val(output_dir)
        val(starting_data)
    output:
    tuple val(sample), path("data/endpoints/${params.project}/$sample/soupX/"), emit: filtered_counts_dir
    script:
    """
    soupX.R \\
    $params.r_lib_path \\
    $sample \\
    $data_type \\
    $params.project \\
    $input_dir \\
    $output_dir \\
    $starting_data
    """
}