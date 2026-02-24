process CONVERT_H5 {
    label 'C2'
    container 'pgc-images.sbgenomics.com/d3b-bixu/swans:v2.1.0'

    input:
    tuple val(meta), path(h5_raw), path(h5_filtered)

    output:
    tuple val(meta), path(sample_id)

    script:
    sample_id = meta.sample_id
    updated_meta = meta + ["input_type": "dir_cellranger"]
    """
    cellranger_h5_to_count_dir.R \\
    --sample $sample_id \\
    --raw $h5_raw \\
    --filtered $h5_filtered
    """
}
