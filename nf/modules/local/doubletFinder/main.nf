process DOUBLETFINDER {
    label 'M6'
    container "francothyroidlab/swans:latest"

    input:
        val(sample)
        val(starting_data)
        path(input_dir)
        val(output_dir)
        val(mito_fraction)
        val(min_feature_threshold)
        val(int_components)
        val(organism)
    output:
    path("data/endpoints/${params.project}/$sample/doubletFinder/tables/*_doublet_ids.txt")
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
    $params.dbl_threads
    """
}