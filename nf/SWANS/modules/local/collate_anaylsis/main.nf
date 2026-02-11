process COLLATE_ANALYSIS{
    label 'C4'
    container "pgc-images.sbgenomics.com/d3b-bixu/swans:v2.1.0"
    input:
        path(input_dirs)
    output:
        path("*")
    script:
    // iterate through samples and dirs to create desired centralized folder structure
    def cp_dirs = ""
    def collate_dir = "data/endpoints/${params.project}/analysis/"
    input_dirs.each { dir ->
        cp_dirs += "cp -r $dir/* $collate_dir;"
    }
    """
    mkdir -p $collate_dir

    $cp_dirs

    cp /SWANS/src/rmd/Interactive_report.Rmd $collate_dir/report/
    """

}