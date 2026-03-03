#!/usr/bin/env nextflow

workflow {
    main:
    user_analyzed_seurat_object = channel.fromPath(params.user_analyzed_seurat_object)
    cluster_annotation_file = channel.fromPath(params.cluster_annotation_file)
    final_user_gene_file = channel.fromPath(params.final_user_gene_file)

    meta = channel.value(
        [
            r_lib_path: params.r_lib_path,
            project: params.project,
            provide_analyzed_seurat_object: "y",
            user_analyzed_seurat_object_meta_sample: params.user_analyzed_seurat_object_meta_sample,
            user_analyzed_seurat_object_meta_experiment: params.user_analyzed_seurat_object_meta_experiment,
            user_analyzed_seurat_object_meta_annotation: params.user_analyzed_seurat_object_meta_annotation,
            user_umap_reduction: params.user_umap_reduction,
            user_tnse_reduction: params.user_tnse_reduction,
            organism: params.organism,
            annotate_provided_final_seurat_object: params.annotate_provided_final_seurat_object,
            min_pct: params.min_pct,
            avg_log2fc_threshold: params.avg_log2fc_threshold,
            final_filtering_threshold: params.final_filtering_threshold,
            final_seurat_normalization_method: params.final_seurat_normalization_method,
            final_seurat_integration_method: params.final_seurat_integration_method,
            final_resolution: params.final_resolution,
            final_storage: params.final_storage,
            final_visualization: params.final_visualization,
            final_conserved_genes: params.final_conserved_genes,
            final_threads: params.final_threads,
            run_trajectory_analysis: params.run_trajectory_analysis,
            partition_trajectory: params.partition_trajectory
        ]
    )
}