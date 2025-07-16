process DOUBLETFINDER {
    label 'M6'
    container "francothyroidlab/swans:latest"

    input:
        val(sample)
        val(starting_data)
        path(input_dir)
        path(output_dir)
        val(mito_fraction)
        val(min_feature_threshold)
        val(int_components)
        val(organism)
    output:
    path('*_doublet_ids.txt')
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
    /usr/local/lib/R/site-library/
    """
}