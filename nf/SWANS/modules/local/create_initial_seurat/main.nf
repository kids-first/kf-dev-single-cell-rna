process CREATE_INITIAL_SEURAT {
    label 'C4'
    container "swans:alpha"

    input:
        tuple val(meta_config), path(collated_data)
        val(sample)
        val(condition)
        path(input_dir)
        val(seurat_creation_source)
        val(run_doubletfinder)
        val(seurat_file_name)
    output:
        path("create_initial_seurat_output"),  emit: analysis_dir
        tuple val(meta_config), path("create_initial_seurat_output/RDS/*.qs"), emit: seurat_file
    script:
    // Create input file table
    def sample_list_str = "samples\tcondition\tpath_to_starting_data\n"
    sample.eachWithIndex { s, i->
        sample_list_str += "${s}\t${condition[i]}\t${input_dir[i]}\n";
        }
    def qc_html_fn = "data/endpoints/$params.project/analysis/report/qc_report/${params.project}_qc_report.html"
    def meta_config_str = ""
    meta_config.each { k, v -> meta_config_str += "${k}: ${v}\n" }

    """
    echo -e "$sample_list_str" > samples.sample_list

    echo -e "$meta_config_str" > prelim_configs.yaml
    
    create_initial_seurat.R \\
    --sample_file samples.sample_list \\
    --project $meta_config.PROJECT \\
    --organism $meta_config.ORGANISM \\
    --seurat_creation_source $seurat_creation_source \\
    --run_doubletfinder $run_doubletfinder \\
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
    params = list(project = "$meta_config.PROJECT", root_dir = "./", data_dir = "./data/endpoints", qc_config = "prelim_configs.yaml"))'

    cp prelim_configs.yaml data/endpoints/$meta_config.PROJECT/analysis/report/

    cp samples.sample_list data/endpoints/$meta_config.PROJECT/analysis/report/${meta_config.PROJECT}_samples.sample_list

    cp -r data/endpoints/${meta_config.PROJECT}/analysis ./create_initial_seurat_output
    """
}