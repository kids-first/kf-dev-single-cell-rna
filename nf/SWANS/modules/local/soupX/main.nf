process SOUPX {
    label 'C4'
    container "pgc-images.sbgenomics.com/d3b-bixu/swans:v2.1.0"

    input:
        val(meta_config)
        tuple val(src), val(sample), path(input_dir)
    output:
    tuple val(sample), path("${sample}_soupX/")

    script:
    def soupX_output_path = "${sample}_soupX"
    """
    soupX.R \\
    --sample $sample \\
    --data_type $meta_config.SOUPX_START \\
    --project $meta_config.PROJECT \\
    --soupX_input_path $input_dir \\
    --soupX_output_path $soupX_output_path \\
    --starter_data $src \\
    $meta_config.RPATH
    """
}