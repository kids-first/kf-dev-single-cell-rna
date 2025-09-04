#!/usr/bin/env nextflow

include { DOUBLETFINDER } from './modules/local/doubletFinder/main.nf'
include { SOUPX } from './modules/local/soupX/main.nf'
include { COLLATE_OUTPUTS } from './modules/local/collate_outputs/main.nf'
include { CREATE_INITIAL_SEURAT } from './modules/local/create_initial_seurat/main.nf'
include { ANALYZE_SEURAT_OBJECT } from './modules/local/analyze_seurat_object/main.nf'
include { CREATE_IMAGES_DGE } from './modules/local/create_images_dge/main.nf'
include { COLLATE_ANALYSIS } from './modules/local/collate_anaylsis/main.nf'


workflow {
    main:
    sample_list = Channel.fromList(params.sample_list)
    condition_list = Channel.fromList(params.condition_list).collect()
    input_dir_list = params.input_dir_list ? Channel.fromPath(params.input_dir_list.class == String ? params.input_dir_list.split(',') as List : params.input_dir_list) : Channel.empty()

    // Create meta dict of common inputs to reduce param passing and to mimic snakemake yaml
    meta = [
        PROJECT: params.project,
        ORGANISM: params.organism,
        RPATH: params.r_lib_path,
        RUN_SOUPX: String.valueOf(!params.disable_soupx),
        SOUPX_START: params.data_type,
        STARTING_DATA: params.starting_data,
        RUN_DOUBLETFINDER: String.valueOf(!params.disable_doubletfinder),
        MITO: params.mito_cutoff,
        RIBO: params.ribo_cutoff,
        MIN_FEATURE_THRESHOLD: params.min_feature_threshold,
        MAX_FEATURE_THRESHOLD: params.max_feature_threshold,
        SEURAT_NORMALIZATION_METHOD: params.normalization_config,
        SEURAT_INTEGRATION_METHOD: params.integration_config,
        RESOLUTION: params.resolution_config,
        COMPONENTS: params.int_components,
        MITO_REGRESSION: params.mito_regression,
        RIBO_REGRESSION: params.ribo_regression,
        CELL_CYCLE_REGRESSION: params.cc_regression,
        NUM_VARIABLE_FEATURES: params.num_var_features,
        SCALE_DATA_FEATURES: params.scale_data_features,
        SPLIT_LAYERS_BY: params.split_layers_by,
        REFERENCE_BASED_INTEGRATION: params.ref_based_integration,
        RUN_AZIMUTH: params.run_azimuth,
        RUN_TRANSFERDATA: params.run_transferdata,
        TSNE: params.include_tsne,
        CONSERVED_GENES: params.conserved_genes,
        VISUALIZATION: params.visualization
    ]

    input_list = sample_list.merge(input_dir_list).map { sample, input_dir -> [sample, input_dir]}
    if (!params.disable_doubletfinder){
        DOUBLETFINDER(
            meta,
            input_list
        )
    }
    if (!params.disable_soupx){
        SOUPX(
            meta,
            input_list
        )
    }
    // collate results so that next step can use them all together
    samples = SOUPX.out.map { it[0] }.concat(DOUBLETFINDER.out.map { it[0] }).collect()
    input_dirs = SOUPX.out.map { it[1] }.collect().combine(DOUBLETFINDER.out.map { it[1] }.collect())
    COLLATE_OUTPUTS(
        meta,
        samples,
        input_dirs
    )
    seurat_filename = "data/endpoints/$params.project/analysis/RDS/${params.project}_initial_seurat_object.qs"
    sample_list_flat = sample_list.collect()
    input_dir_list_flat = input_dir_list.collect()
    CREATE_INITIAL_SEURAT(
        COLLATE_OUTPUTS.out,
        sample_list_flat,
        condition_list,
        input_dir_list_flat,
        seurat_filename
    )
    ANALYZE_SEURAT_OBJECT(
        CREATE_INITIAL_SEURAT.out.seurat_file,
        params.aso_memory
    )
    CREATE_IMAGES_DGE(
        params.storage,
        ANALYZE_SEURAT_OBJECT.out.analyzed_seurat_object_file
    )
    analysis_dirs = CREATE_INITIAL_SEURAT.out.analysis_dir.combine(ANALYZE_SEURAT_OBJECT.out.analysis_path).combine(CREATE_IMAGES_DGE.out)
    COLLATE_ANALYSIS(
        analysis_dirs
    )
}