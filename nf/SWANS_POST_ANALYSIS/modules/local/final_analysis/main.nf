process FINAL_ANALYSIS {
    label 'M8'
    container "pgc-images.sbgenomics.com/d3b-bixu/swans:v2.1.0"

    input:
        val(meta_config)
        path(existing_analysis_tar)
        path(cluster_annotation_file)
        path(final_user_gene_file)
    output:
        path("data/endpoints"), emit: analysis_path
    script:
    def rds_root = "data/endpoints/${meta_config.project}/analysis/RDS/${meta_config.project}_analyzed_seurat_object"
    """
    tar -xf $existing_analysis_tar \\
    && mkdir -p data/endpoints/${meta_config.project} && mv analysis data/endpoints/${meta_config.project}/

    if [ -f "${rds_root}.RDS" ];then
        mv ${rds_root}.RDS ${rds_root}.rds
    fi

    final_analysis.R \\
    $meta_config.project \\
    $meta_config.r_lib_path \\
    ${rds_root}.rds \\
    $meta_config.final_seurat_normalization_method \\
    $meta_config.final_seurat_integration_method \\
    $meta_config.final_resolution \\
    $cluster_annotation_file \\
    $meta_config.user_tnse_reduction \\
    $final_user_gene_file \\
    $meta_config.final_threads \\
    $meta_config.final_filtering_threshold \\
    $meta_config.avg2fc \\
    $meta_config.min_pct \\
    $meta_config.final_conserved_genes \\
    $meta_config.organism \\
    $meta_config.final_storage \\
    $meta_config.provide_analyzed_seurat_object \\
    $meta_config.annotate_provided_final_seurat_object \\
    $meta_config.final_visualization \\
    $meta_config.user_analyzed_seurat_object_meta_sample \\
    $meta_config.user_analyzed_seurat_object_meta_experiment \\
    $meta_config.user_analyzed_seurat_object_meta_annotation \\
    $meta_config.user_umap_reduction \\
    $meta_config.user_tnse_reduction \\
    $meta_config.memory_mb

    """

}