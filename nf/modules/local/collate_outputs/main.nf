process COLLATE_OUTPUTS{
    label 'C4'
    container "ubuntu:latest"
    input:
        val(samples)
        path(input_dirs)
        val(root_dir)
    output:
        path("collated/*")
    script:
    // iterate through samples and dirs to create desired centralized folder structure
    def createdirs = ""
    samples.eachWithIndex { sample, index ->
        def cur_dir = input_dirs[index]
        def tool = cur_dir.name.replaceFirst("${sample}_", "")
        def collate_dir = "collated/$root_dir/$params.project/$sample/$tool/"
        createdirs += "mkdir -p $collate_dir && cp -r $cur_dir/* $collate_dir;"
    }
    """
    $createdirs
    """

}