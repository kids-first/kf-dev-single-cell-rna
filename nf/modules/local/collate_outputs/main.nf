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
    def i=0
    samples.each {
        def sample_name = it
        def cur_dir = input_dirs[i]
        def tool = cur_dir.name.replaceFirst("${sample_name}_", "")
        def collate_dir = "collated/$root_dir/$params.project/$sample_name/$tool/"
        createdirs += "mkdir -p $collate_dir && cp -r $cur_dir/* $collate_dir;"
        i += 1;
    }
    """
    $createdirs
    """

}