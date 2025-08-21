process CREATE_IMAGES_DGE {
    label 'M4'
    container "swans:alpha"

    input:
        val(storage)
        val(normalization_method)
        val(integration_method)
        val(resolution)
        val (conserved_genes)
        path(analyzed_seurat_object)
        val(tsne_plot)
        val(report_table_path)
        val(visualization)
    output:
        path("./create_images_dge_output")
    script:
    def args = task.ext.args ?: ''
    """
    create_images_DGE.R \\
    --project $params.project \\
    --storage $storage \\
    --normalization_method $normalization_method \\
    --integration_method $integration_method \\
    --resolution $resolution \\
    --conserved_genes $conserved_genes \\
    --analyzed_seurat_object $analyzed_seurat_object \\
    --processes 4 \\
    --tsne_plot $tsne_plot \\
    --report_table_path $report_table_path \\
    --visualization $visualization \\
    $args \\
    $params.r_lib_path

    cp -r data/endpoints/${params.project}/analysis ./create_images_dge_output
    """
}