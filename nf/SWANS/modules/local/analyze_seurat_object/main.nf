process ANALYZE_SEURAT_OBJECT {
    label 'ASO'
    container "swans:alpha"

    input:
        tuple val(meta_config), path(initial_seurat_object)
        val(aso_memory)
    output:
        path("analyze_seurat_object_output"), emit: analysis_path
        tuple val(meta_config), path("analyze_seurat_object_output/RDS/${meta_config.PROJECT}_analyzed_seurat_object.qs"), emit: analyzed_seurat_object_file
    script:
    def analyzed_seurat_object = "data/endpoints/$meta_config.PROJECT/analysis/RDS/${meta_config.PROJECT}_analyzed_seurat_object.qs"
    def report_path_figures = "data/endpoints/$params.project/analysis/report/figures" 
    def args = task.ext.args ?: ''
    """
    echo -e "$aso_memory" > memory.txt
    
    analyze_seurat_object.R \\
    --initial_seurat_object $initial_seurat_object \\
    --project $meta_config.PROJECT \\
    --organism $meta_config.ORGANISM \\
    --n_components $meta_config.COMPONENTS \\
    --mito_regression $meta_config.MITO_REGRESSION \\
    --ribo_regression $meta_config.RIBO_REGRESSION \\
    --cc_regression $meta_config.CELL_CYCLE_REGRESSION \\
    --num_var_features $meta_config.NUM_VARIABLE_FEATURES \\
    --scale_data_features $meta_config.SCALE_DATA_FEATURES \\
    --split_layers_by $meta_config.SPLIT_LAYERS_BY \\
    --normalization_config $meta_config.SEURAT_NORMALIZATION_METHOD \\
    --integration_config $meta_config.SEURAT_INTEGRATION_METHOD \\
    --ref_based_integration $meta_config.REFERENCE_BASED_INTEGRATION \\
    --run_azimuth $meta_config.RUN_AZIMUTH \\
    --run_transferdata $meta_config.RUN_TRANSFERDATA \\
    --resolution_config $meta_config.RESOLUTION \\
    --include_tsne $meta_config.TSNE \\
    --analyzed_seurat_object $analyzed_seurat_object \\
    --report_path_figures $report_path_figures \\
    --processes $task.cpus \\
    --memory_file memory.txt \\
    $args \\
    $meta_config.RPATH

    cp -r data/endpoints/${meta_config.PROJECT}/analysis ./analyze_seurat_object_output
    """
}