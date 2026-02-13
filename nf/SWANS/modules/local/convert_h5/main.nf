process CONVERT_H5 {
    label 'C2'
    container 'pgc-images.sbgenomics.com/d3b-bixu/swans:v2.1.0'

    input:
    tuple val(sample), val(src), path(h5_files)

    output:
    path(sample)

    script:
    
    def raw_h5_file = src[0] == "h5_raw" ? h5_files[0] : h5_files[1]
    def filtered_h5_file = src[0] == "h5_filtered" ? h5_files[0] : h5_files[1]
    """
    cellranger_h5_to_count_dir.R \\
    --sample $sample \\
    --raw $raw_h5_file \\
    --filtered $filtered_h5_file
    """
}