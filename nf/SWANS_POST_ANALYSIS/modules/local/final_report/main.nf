process FINAL_REPORT {
    label 'C8'
    container "pgc-images.sbgenomics.com/d3b-bixu/swans:v2.1.0"

    input:
        val(meta_config)
        path(sample_list)
        path(prelim_config)
        path(final_analysis_dir)
        path(cluster_annotation_file)
        path(final_user_gene_file)

    output:
        path("${meta_config.PROJECT}_final_report.html")
    script:
    def script = "/SWANS/src/rmd/final_report.Rmd"
    def out_html = "${meta_config.PROJECT}_final_report.html"
    def meta_config_str = meta_config.collect{ k, v -> "${k}: ${v}" }.join('\n')
    """
    echo -e "$meta_config_str" > post_annotation_config.yaml

    echo -e "CLUSTER_ANNOTATION_FILE: $cluster_annotation_file" >> post_annotation_config.yaml

    echo -e "FINAL_USER_GENE_FILE: $final_user_gene_file" >> post_annotation_config.yaml

    cp $script ./

    Rscript -e 'library(rmarkdown); \
    rmarkdown::render("final_report.Rmd", \
    output_file="$out_html", \
    params = list(root_dir = "./", data_dir = "./data/endpoints/", prelim_qc_config = "prelim_config.yaml", post_annotation_config = "post_annotation_config.yaml")
    )'

    """

}