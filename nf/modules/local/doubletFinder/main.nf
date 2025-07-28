process DOUBLETFINDER {
    label 'M6'
    container "swans:alpha"

    input:
        val(sample)
        path(input_dir)
        val(starting_data)
        val(output_dir)
        val(mito_fraction)
        val(min_feature_threshold)
        val(int_components)
        val(organism)
    output:
    path("${sample}_doubletFinder/")
    script:
    """
    doubletFinder.R \\
    $sample \\
    $params.project \\
    $starting_data \\
    $input_dir \\
    $output_dir \\
    $mito_fraction \\
    $min_feature_threshold \\
    $int_components \\
    $organism \\
    $params.r_lib_path \\
    $params.dbl_threads \\
    && mkdir ${sample}_doubletFinder \\
    && mv $output_dir/$params.project/$sample/* ${sample}_doubletFinder/
    """
}