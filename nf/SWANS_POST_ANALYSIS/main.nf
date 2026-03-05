#!/usr/bin/env nextflow

include { FINAL_ANALYSIS } from './modules/local/final_analysis/main.nf'
include { TRAJECTORY_ANALYSIS } from './modules/local/trajectory_analysis/main.nf'
include { FINAL_REPORT } from './modules/local/final_report/main.nf'
include { TAR_OUTPUTS } from './modules/local/tar/main.nf'

workflow {
    main:
    existing_analysis_tar = channel.fromPath(params.existing_analysis_tar)
    sample_list = channel.fromPath(params.sample_list)
    prelim_config = channel.fromPath(params.prelim_config)
    cluster_annotation_file = channel.fromPath(params.cluster_annotation_file)
    final_user_gene_file = channel.fromPath(params.final_user_gene_file)

    meta = channel.value(
        [
            RPATH: params.r_lib_path,
            PROJECT: params.project,
            PROVIDE_ANALYZED_SEURAT_OBJECT: "y",
            USER_ANALYZED_SEURAT_OBJECT_META_SAMPLE: params.user_analyzed_seurat_object_meta_sample,
            USER_ANALYZED_SEURAT_OBJECT_META_EXPERIMENT: params.user_analyzed_seurat_object_meta_experiment,
            USER_ANALYZED_SEURAT_OBJECT_META_ANNOTATION: params.user_analyzed_seurat_object_meta_annotation,
            USER_UMAP_REDUCTION: params.user_umap_reduction,
            USER_TNSE_REDUCTION: params.user_tnse_reduction,
            ORGANISM: params.organism,
            RUN_TRAJECTORY_ANALYSIS: params.run_trajectory_analysis ? 'y' : 'n',
            ANNOTATE_PROVIDED_FINAL_SEURAT_OBJECT: params.annotate_provided_final_seurat_object,
            MIN_PCT: params.min_pct,
            AVG_LOG2FC_THRESHOLD: params.avg_log2fc_threshold,
            FINAL_FILTERING_THRESHOLD: params.final_filtering_threshold,
            FINAL_SEURAT_NORMALIZATION_METHOD: params.final_seurat_normalization_method,
            FINAL_SEURAT_INTEGRATION_METHOD: params.final_seurat_integration_method,
            FINAL_RESOLUTION: params.final_resolution,
            FINAL_STORAGE: params.final_storage,
            FINAL_VISUALIZATION: params.final_visualization,
            FINAL_CONSERVED_GENES: params.final_conserved_genes,
            FINAL_THREADS: params.final_threads,
            PARTITION_TRAJECTORY: params.partition_trajectory,
            MEMORY_MB: params.memory_mb
        ]
    )

    FINAL_ANALYSIS(
        meta,
        existing_analysis_tar,
        cluster_annotation_file,
        final_user_gene_file
    )
    if (params.run_trajectory_analysis){
        TRAJECTORY_ANALYSIS(
            meta,
            FINAL_ANALYSIS.out.analysis_path,
            cluster_annotation_file
        )
    }
    final_report_input_dir = TRAJECTORY_ANALYSIS.out ?: FINAL_ANALYSIS.out.analysis_path
    FINAL_REPORT(
        meta,
        sample_list,
        prelim_config,
        final_report_input_dir,
        cluster_annotation_file,
        final_user_gene_file
    )

    tar_output_input = final_report_input_dir.map { dir_path -> ["data/endpoints/$params.project", "${params.project}.final_analysis.tar.gz", dir_path] }
    TAR_OUTPUTS(
        tar_output_input
    )

}