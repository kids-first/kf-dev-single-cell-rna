process COLLATE_OUTPUTS{
    label 'C4'
    container "ubuntu:latest"
    input:
        val(meta_config)
        val(meta)
        path(input_dirs)
    output:
        tuple val(meta_config), path("collated/*")
    script:
    // iterate through samples and dirs to create desired centralized folder structure
    def createdirs = ""
    meta.eachWithIndex { metadata, index ->
        def sample = metadata.containsKey("remap") ? metadata.remap : metadata.sample_id
        def cur_dir = input_dirs[index]
        def tool = metadata.input_type.contains("cellranger") ? "cellranger" : metadata.input_type.replaceAll("tar_|dir_", "")
        def collate_dir = "collated/data/endpoints/$params.project/$sample/$tool/"
        createdirs += "mkdir -p $collate_dir && cp -r $cur_dir/* $collate_dir;"
    }
    """
    $createdirs
    """

}