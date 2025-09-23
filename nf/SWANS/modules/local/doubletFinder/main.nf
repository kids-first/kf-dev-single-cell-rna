process DOUBLETFINDER {
    label 'M6'
    container "swans:alpha"

    input:
        val(meta_config)
        tuple val(sample), path(input_dir)

    output:
    tuple val(sample), path("${sample}_doubletFinder/")
    script:
    def doubletFinder_output_path = "${sample}_doubletFinder/"
    """
    doubletFinder.R \\
    --sample $sample \\
    --project $meta_config.PROJECT \\
    --starting_data $meta_config.STARTING_DATA \\
    --input_path $input_dir \\
    --output_path $doubletFinder_output_path \\
    --mito_cutoff $meta_config.MITO \\
    --min_feature_threshold $meta_config.MIN_FEATURE_THRESHOLD \\
    --components $meta_config.COMPONENTS \\
    --organism $meta_config.ORGANISM \\
    --processes $task.cpus \\
    $meta_config.RPATH
    """
}