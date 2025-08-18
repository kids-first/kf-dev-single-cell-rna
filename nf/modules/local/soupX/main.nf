process SOUPX {
    label 'C4'
    container "swans:alpha"

    input:
        tuple val(sample), path(input_dir)
        val(data_type)
        val(output_dir)
        val(starting_data)
    output:
    tuple val(sample), path("${sample}_soupX/")

    script:
    def soupX_output_path = "$output_dir/$params.project/$sample/soupX/"
    """
    soupX.R \\
    --sample $sample \\
    --data_type $data_type \\
    --project $params.project \\
    --soupX_input_path $input_dir \\
    --soupX_output_path $soupX_output_path \\
    --starter_data $starting_data \\
    $params.r_lib_path \\
    && mkdir ${sample}_soupX \\
    && mv $soupX_output_path/* ${sample}_soupX/
    """
}