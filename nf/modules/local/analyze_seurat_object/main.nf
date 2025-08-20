process ANALYZE_SEURAT_OBJECT {
    label 'ASO'
    container "swans:alpha"

    input:
        path(initial_seurat_object)
        val(n_components)
        val(mito_regression)
        val(ribo_regression)
        val(cc_regression)
        val(num_var_features)
        val(scale_data_features)
        val(split_layers_by)
        val(normalization_config)
        val(integration_config)
        val(ref_based_integration)
        val(run_azimuth)
        val(run_transferdata)
        val(resolution_config)
        val(include_tsne)
        val(analyzed_seurat_object)
        val(report_path_figures)
    output:
        path("data/endpoints/${params.project}/analysis"), emit: anaylsis_path
        path(analyzed_seurat_object), emit: analyzed_seurat_object_file
    script:
    def args = task.ext.args ?: ''
    """
    echo -e "$params.aso_memory" > memory.txt
    
    analyze_seurat_object.R \\
    --initial_seurat_object $initial_seurat_object \\
    --project $params.project \\
    --organism $params.organism \\
    --n_components $n_components \\
    --mito_regression $mito_regression \\
    --ribo_regression $ribo_regression \\
    --cc_regression $cc_regression \\
    --num_var_features $num_var_features \\
    --scale_data_features $scale_data_features \\
    --split_layers_by $split_layers_by \\
    --normalization_config $normalization_config \\
    --integration_config $integration_config \\
    --ref_based_integration $ref_based_integration \\
    --run_azimuth $run_azimuth \\
    --run_transferdata $run_transferdata \\
    --resolution_config $resolution_config \\
    --include_tsne $include_tsne \\
    --analyzed_seurat_object $analyzed_seurat_object \\
    --report_path_figures $report_path_figures \\
    --processes $params.aso_processes \\
    --memory_file memory.txt \\
    $args \\
    $params.r_lib_path

    cp -r data/endpoints/${params.project}/analysis ./
    """
}