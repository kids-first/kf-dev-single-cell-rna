process CREATE_INITIAL_SEURAT {
    label 'C4'
    container "pgc-images.sbgenomics.com/d3b-bixu/swans:v2.1.0"

    input:
        tuple val(meta_config), path(collated_data)
        val(sample)
        val(condition)
        path(input_dir)
        val(seurat_file_name)
    output:
        tuple val(meta_config), path("${meta_config.PROJECT}_initial_seurat_object.qs"), emit: seurat_qs
        tuple val(meta_config), path("${meta_config.PROJECT}_qc_report.html"), emit: qc_report
    script:
    // Create input file table
    def sample_list_str = "samples\tcondition\tpath_to_starting_data\n"
    sample.eachWithIndex { s, i->
        sample_list_str += "${s}\t${condition[i]}\t${input_dir[i]}\n";
        }
    def qc_html_fn = "${params.project}_qc_report.html"
    def meta_config_str = ""
    meta_config.each { k, v -> meta_config_str += "${k}: ${v}\n" }

    """
    echo -e "$sample_list_str" > samples.sample_list

    echo -e "$meta_config_str" > prelim_configs.yaml
    
    create_initial_seurat.R \\
    --sample_file samples.sample_list \\
    --project $meta_config.PROJECT \\
    --organism $meta_config.ORGANISM \\
    --seurat_creation_source ${meta_config.RUN_SOUPX == "true" ? "soupX": "cellranger" } \\
    --run_doubletfinder ${meta_config.RUN_DOUBLETFINDER == "true" ? "y" : "n"} \\
    --mito_cutoff $meta_config.MITO \\
    --ribo_cutoff $meta_config.RIBO \\
    --min_feature_threshold $meta_config.MIN_FEATURE_THRESHOLD \\
    --max_feature_threshold $meta_config.MAX_FEATURE_THRESHOLD \\
    --seurat_file_name $seurat_file_name \\
    $meta_config.RPATH
    
    echo "Generating R markdown QC report"
    
    cp /SWANS/src/rmd/qc_report.Rmd ./
    
    Rscript -e 'library(rmarkdown); \
    rmarkdown::render("qc_report.Rmd", \
    output_file="$qc_html_fn", \
    params = list(root_dir = "./", data_dir = "./data/endpoints", qc_config = "prelim_configs.yaml"))'

    """
}