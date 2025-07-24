process CREATE_INITIAL_SEURAT {
    label 'C4'
    container "swans:alpha"

    input:
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
    output:
    tuple val(sample), path("data/endpoints/${params.project}/$sample/doubletFinder/tables/*_doublet_ids.txt")
    script:
    // Create input file table
    def sample_list_str = ""
    def i = 0
    sample.each { s->
        sample_list_str += s + "\t" + condition[i] + "\t" + input_dir[i] + "\n";
        i += 1
        }
    """
    echo -e "$sample_list_str" > samples.sample_list
    create_initial_seurat.R \\
    $sample \\
    $params.project \\
    $seurat_creation_source \\
    $input_dir \\
    $run_doubletfinder \\
    $mito_cutoff \\
    $ribo_cutoff \\
    $min_feature_threshold \\
    $max_feature_threshold
    $seurat_file_name \\
    $params.r_lib_path
    """
}