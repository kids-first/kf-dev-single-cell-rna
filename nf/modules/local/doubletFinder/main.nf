process DOUBLETFINDER {
    label 'M6'
    container "swans:alpha"

    input:
        tuple val(sample), path(input_dir)
        val(starting_data)
        val(output_dir)
        val(mito_fraction)
        val(min_feature_threshold)
        val(int_components)
        val(organism)
    output:
    tuple val(sample), path("${sample}_doubletFinder/")
    script:
    """
    doubletFinder.R \\
    --sample $sample \\
    --project $params.project \\
    --starting_data $starting_data \\
    --input_path $input_dir \\
    --output_path $output_dir \\
    --mito_cutoff $mito_fraction \\
    --min_feature_threshold $min_feature_threshold \\
    --components $int_components \\
    --organism $organism \\
    --processes $task.cpus \\
    $params.r_lib_path \\
    && mkdir ${sample}_doubletFinder \\
    && mv $output_dir/$params.project/$sample/doubletFinder/* ${sample}_doubletFinder/
    """
}