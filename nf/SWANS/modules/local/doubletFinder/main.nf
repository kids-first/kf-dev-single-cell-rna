process DOUBLETFINDER {
    label 'M6'
    container "pgc-images.sbgenomics.com/d3b-bixu/swans:v2.1.0"

    input:
        val(meta_config)
        tuple val(meta), path(input_dir)

    output:
    tuple val(meta), path("${meta.sample_id}_doubletFinder/")
    script:
    def src = meta.input_type.contains("cellranger") ? "cellranger" : "matrix"
    def doubletFinder_output_path = "${meta.sample_id}_doubletFinder/"
    """
    doubletFinder.R \\
    --sample $meta.sample_id \\
    --project $meta_config.PROJECT \\
    --starting_data $src \\
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