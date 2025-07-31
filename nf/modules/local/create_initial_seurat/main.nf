process CREATE_INITIAL_SEURAT {
    label 'C4'
    container "swans:alpha"

    input:
        path(collated_data)
        val(sample)
        val(condition)
        path(input_dir)
        val(seurat_creation_source)
        val(run_doubletfinder)
        val(mito_cutoff)
        val(ribo_cutoff)
        val(min_feature_threshold)
        val(max_feature_threshold)
        val(seurat_file_name)
        path(qc_config)
    output:
    path("data/endpoints/${params.project}/analysis")
    script:
    // Create input file table
    def sample_list_str = "samples\tcondition\tpath_to_starting_data\n"
    def i = 0
    sample.each { s->
        sample_list_str += s + "\t" + condition[i] + "\t" + input_dir[i] + "\n";
        i += 1
        }
    def qc_html_fn = "data/endpoints/$params.project/analysis/report/qc_report/${params.project}_qc_report.html"
    """
    echo -e "$sample_list_str" > samples.sample_list
    create_initial_seurat.R \\
    --sample_file samples.sample_list \\
    --project $params.project \\
    --organism $params.organism \\
    --seurat_creation_source $seurat_creation_source \\
    --run_doubletfinder $run_doubletfinder \\
    --mito_cutoff $mito_cutoff \\
    --ribo_cutoff $ribo_cutoff \\
    --min_feature_threshold $min_feature_threshold \\
    --max_feature_threshold $max_feature_threshold \\
    --seurat_file_name $seurat_file_name \\
    $params.r_lib_path \\
    && echo "Generating R markdown QC report" \\
    && cp /SWANS/src/rmd/qc_report.Rmd ./ \\
    && Rscript -e 'library(rmarkdown); rmarkdown::render("qc_report.Rmd", output_file=\"$qc_html_fn\", params = list(root_dir = "./", data_dir = "./data/endpoints/", qc_config = "$qc_config.name"))' \\
    """
}