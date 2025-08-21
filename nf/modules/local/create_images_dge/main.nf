process CREATE_IMAGES_DGE {
    label 'M4'
    container "swans:alpha"

    input:
        val(storage)
        tuple val(meta_config), path(analyzed_seurat_object)
    output:
        path("./create_images_dge_output")
    script:
    def report_table_path = "data/endpoints/$meta_config.PROJECT/analysis/final_analysis/tables"
    def args = task.ext.args ?: ''
    """
    create_images_DGE.R \\
    --project $meta_config.PROJECT \\
    --storage $storage \\
    --normalization_method $meta_config.SEURAT_NORMALIZATION_METHOD \\
    --integration_method $meta_config.SEURAT_INTEGRATION_METHOD \\
    --resolution $meta_config.RESOLUTION \\
    --conserved_genes $meta_config.CONSERVED_GENES \\
    --analyzed_seurat_object $analyzed_seurat_object \\
    --processes $task.cpus \\
    --tsne_plot $meta_config.TSNE \\
    --report_table_path $report_table_path \\
    --visualization $meta_config.VISUALIZATION \\
    $args \\
    $meta_config.RPATH

    cp -r data/endpoints/$meta_config.PROJECT/analysis ./create_images_dge_output
    """
}