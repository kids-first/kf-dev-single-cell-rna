process TRAJECTORY_ANALYSIS {
    label 'C4'
    container "pgc-images.sbgenomics.com/d3b-bixu/swans:v2.1.0"

    input:
        val(meta_config)
        path(final_analysis_dir)
        path(cluster_annotation_file)
    output:
        path("data")
    script:
    def final_rds_loc = "data/endpoints/${meta_config.PROJECT}/analysis/RDS/${meta_config.PROJECT}_final_analyzed_seurat_object.qs"
    """

    trajectory_analysis.R \\
    $meta_config.PROJECT \\
    $final_rds_loc \\
    $cluster_annotation_file \\
    $meta_config.FINAL_SEURAT_NORMALIZATION_METHOD \\
    $meta_config.FINAL_SEURAT_INTEGRATION_METHOD \\
    $meta_config.FINAL_RESOLUTION \\
    $meta_config.PARTITION_TRAJECTORY \\
    $meta_config.FINAL_STORAGE \\
    $meta_config.RPATH \\
    $meta_config.PROVIDE_ANALYZED_SEURAT_OBJECT \\
    $meta_config.ANNOTATE_PROVIDED_FINAL_SEURAT_OBJECT

    """

}