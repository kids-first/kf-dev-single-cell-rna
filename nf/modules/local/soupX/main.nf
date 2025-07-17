process SOUPX {
    label 'M4'
    container "francothyroidlab/swans:latest"

    input:
        val(sample)
        val(data_type)
        path(input_dir)
        val(output_dir)
    output:
    path("data/endpoints/${params.project}/$sample/soupX/barcodes.tsv.gz"), emit: barcodes_path
    path("data/endpoints/${params.project}/$sample/soupX/features.tsv.gz"), emit: features_path
    path("data/endpoints/${params.project}/$sample/soupX/matrix.tsv.gz"), emit: matrix_path
    script:
    """
    soupX.R \\
    $params.r_lib_path \\
    $sample \\
    $data_type \\
    $params.project \\
    $input_dir \\
    $output_dir
    """
}