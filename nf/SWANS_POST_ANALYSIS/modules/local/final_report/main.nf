process FINAL_REPORT {
    label 'C8'
    container "pgc-images.sbgenomics.com/d3b-bixu/swans:v2.1.0"

    input:
        val(meta_config)
        path(final_analysis_dir)
    output:
        path("data")
    script:
    def script = "/SWANS/src/rmd/final_report.Rmd"
    def out_html = "data/endpoints/${meta_config.project}/analysis/final_analysis/${meta_config.project}_final_report.html"
    """
    complete_path=\$PWD/$out_html

    echo \$complete_path

    Rscript -e 'library(rmarkdown); \
    rmarkdown::render("${script}",\
    output_file="\$(complete_path)")'

    """

}