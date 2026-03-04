process TRAJECTORY_ANALYSIS {
    label 'C8'
    container "pgc-images.sbgenomics.com/d3b-bixu/swans:v2.1.0"

    input:
        val(meta_config)
        path(final_analysis_dir)
        path(cluster_annotation_file)
    output:
        path("data")
    script:
    def final_rds_loc = "data/endpoints/${meta_config.project}/analysis/RDS/${meta_config.project}_final_analyzed_seurat_object.rds"
    """

    trajectory_analysis.R \\
    $meta_config.project \\
    $final_rds_loc \\
    $cluster_annotation_file \\
    $meta_config.final_seurat_normalization_method \\
    $meta_config.final_seurat_integration_method \\
    $meta_config.final_resolution \\
    $meta_config.partition_trajectory \\
    $meta_config.final_storage \\
    $meta_config.r_lib_path \\
    $meta_config.provide_analyzed_seurat_object \\
    $meta_config.annotate_provided_final_seurat_object

    """

}