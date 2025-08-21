#!/usr/bin/env nextflow

include { DOUBLETFINDER } from './modules/local/doubletFinder/main.nf'
include { SOUPX } from './modules/local/soupX/main.nf'
include { COLLATE_OUTPUTS } from './modules/local/collate_outputs/main.nf'
include { CREATE_INITIAL_SEURAT } from './modules/local/create_initial_seurat/main.nf'
include  {ANALYZE_SEURAT_OBJECT } from './modules/local/analyze_seurat_object/main.nf'
include { CREATE_IMAGES_DGE } from './modules/local/create_images_dge/main.nf'
include { COLLATE_ANALYSIS } from './modules/local/collate_anaylsis/main.nf'


def create_yml_config() {
    def yml_config_str = """\
        PROJECT: ${params.project}
        ORGANISM: ${params.organism}
        RUN_SOUPX: ${String.valueOf(!params.disable_soupx)}
        RUN_DOUBLETFINDER: ${String.valueOf(!params.disable_doubletfinder)}
        MITO: ${params.mito_cutoff}
        RIBO: ${params.ribo_cutoff}
        MIN_FEATURE_THRESHOLD: ${params.min_feature_threshold}
        MAX_FEATURE_THRESHOLD: ${params.max_feature_threshold}\
    """.stripIndent()
    def yml_config = file('qc_config.yml')
    yml_config.text = yml_config_str
    return yml_config
}

workflow {
    main:
    sample_list = Channel.fromList(params.sample_list)
    condition_list = Channel.fromList(params.condition_list).collect()
    starting_data = Channel.value(params.starting_data)
    input_dir_list = Channel.fromPath(params.input_dir_list, type: 'dir')
    data_dir = params.data_dir

    input_list = sample_list.merge(input_dir_list).map { sample, input_dir -> [sample, input_dir]}
    if (!params.disable_doubletfinder){
        DOUBLETFINDER(
            input_list,
            starting_data,
            data_dir,
            params.mito_cutoff,
            params.min_feature_threshold,
            params.int_components,
            params.organism
        )
    }
    if (!params.disable_soupx){
        SOUPX(
            input_list,
            params.data_type,
            data_dir,
            starting_data
        )
    }
    // collate results so that next step can use them all together
    samples = SOUPX.out.map { it[0] }.concat(DOUBLETFINDER.out.map { it[0] }).collect()
    input_dirs = SOUPX.out.map { it[1] }.collect().combine(DOUBLETFINDER.out.map { it[1] }.collect())
    COLLATE_OUTPUTS(
        samples,
        input_dirs,
        data_dir
    )
    COLLATE_OUTPUTS.out.view()
    seurat_filename = "data/endpoints/$params.project/analysis/RDS/${params.project}_initial_seurat_object.qs"
    qc_config = create_yml_config()
    print qc_config.text
    CREATE_INITIAL_SEURAT(
        COLLATE_OUTPUTS.out,
        sample_list.collect(),
        condition_list,
        input_dir_list.collect(),
        params.disable_soupx ? "cellranger" : "soupX",
        'y',
        params.mito_cutoff,
        params.ribo_cutoff,
        params.min_feature_threshold,
        params.max_feature_threshold,
        seurat_filename,
        qc_config
    )
    ANALYZE_SEURAT_OBJECT(
        CREATE_INITIAL_SEURAT.out.seurat_file,
        params.int_components,
        params.mito_regression,
        params.ribo_regression,
        params.cc_regression,
        params.num_var_features,
        params.scale_data_features,
        params.split_layers_by,
        params.normalization_config,
        params.integration_config,
        params.ref_based_integration,
        params.run_azimuth,
        params.run_transferdata,
        params.resolution_config,
        params.include_tsne,
        params.analyzed_seurat_object,
        params.report_path_figures
    )
    report_table_path = "${params.data_dir}/${params.project}/analysis/final_analysis/tables"
    CREATE_IMAGES_DGE(
        params.storage,
        params.normalization_config,
        params.integration_config,
        params.resolution_config,
        params.conserved_genes,
        ANALYZE_SEURAT_OBJECT.out.analyzed_seurat_object_file,
        params.include_tsne,
        report_table_path,
        params.visualization
    )
    analysis_dirs = CREATE_INITIAL_SEURAT.out.analysis_dir.combine(ANALYZE_SEURAT_OBJECT.out.analysis_path).combine(CREATE_IMAGES_DGE.out)
    COLLATE_ANALYSIS(
        analysis_dirs
    )
}