process SOUPX {
    label 'C4'
    container "pgc-images.sbgenomics.com/d3b-bixu/swans:v2.1.0"

    input:
        val(meta_config)
        tuple val(meta), path(input_dir)
    output:
    tuple val(meta), path("${meta.sample_id}_soupX/")

    script:
    def soupX_output_path = "${meta.sample_id}_soupX"
    """
    soupX.R \\
    --sample $meta.sample_id \\
    --data_type $meta_config.SOUPX_START \\
    --project $meta_config.PROJECT \\
    --soupX_input_path $input_dir \\
    --soupX_output_path $soupX_output_path \\
    --starter_data cellranger \\
    $meta_config.RPATH
    """
}