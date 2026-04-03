process FINAL_ANALYSIS {
    label 'M8'
    container "pgc-images.sbgenomics.com/d3b-bixu/swans:v2.1.0"

    input:
        val(meta_config)
        path(existing_analysis_tar)
        path(cluster_annotation_file)
        path(final_user_gene_file)
    output:
        path("data"), emit: analysis_path
        path("data/endpoints/${meta_config.PROJECT}/analysis/report/prelim_configs.yaml"), emit: extracted_prelim_config
    script:
    def rds_root = "data/endpoints/${meta_config.PROJECT}/analysis/RDS/${meta_config.PROJECT}_analyzed_seurat_object"
    """
    tar xf $existing_analysis_tar \\
    && mkdir -p data/endpoints/${meta_config.PROJECT} && mv analysis data/endpoints/${meta_config.PROJECT}/

    final_analysis.R \\
    $meta_config.PROJECT \\
    $meta_config.RPATH \\
    ${rds_root}.qs \\
    $meta_config.FINAL_SEURAT_NORMALIZATION_METHOD \\
    $meta_config.FINAL_SEURAT_INTEGRATION_METHOD \\
    $meta_config.FINAL_RESOLUTION \\
    $cluster_annotation_file \\
    $meta_config.USER_TNSE_REDUCTION \\
    $final_user_gene_file \\
    $meta_config.FINAL_THREADS \\
    $meta_config.FINAL_FILTERING_THRESHOLD \\
    $meta_config.AVG_LOG2FC_THRESHOLD \\
    $meta_config.MIN_PCT \\
    $meta_config.FINAL_CONSERVED_GENES \\
    $meta_config.ORGANISM \\
    $meta_config.FINAL_STORAGE \\
    $meta_config.PROVIDE_ANALYZED_SEURAT_OBJECT \\
    $meta_config.ANNOTATE_PROVIDED_FINAL_SEURAT_OBJECT \\
    $meta_config.FINAL_VISUALIZATION \\
    $meta_config.USER_ANALYZED_SEURAT_OBJECT_META_SAMPLE \\
    $meta_config.USER_ANALYZED_SEURAT_OBJECT_META_EXPERIMENT \\
    $meta_config.USER_ANALYZED_SEURAT_OBJECT_META_ANNOTATION \\
    $meta_config.USER_UMAP_REDUCTION \\
    $meta_config.USER_TNSE_REDUCTION \\
    $meta_config.MEMORY_MB

    """

}