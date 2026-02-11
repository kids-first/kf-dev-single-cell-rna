process CONVERT_H5 {
    label 'C2'
    container 'pgc-images.sbgenomics.com/d3b-bixu/swans:v2.1.0'

    input:
    tuple val(sample), path(raw_h5_file), path(filtered_h5_file)

    output:
    path(sample)

    script:
    """
    cellranger_h5_to_count_dir.R \\
    --sample $sample \\
    --raw_h5 $raw_h5_file \\
    --filtered_h5 $filtered_h5_file
    """
}