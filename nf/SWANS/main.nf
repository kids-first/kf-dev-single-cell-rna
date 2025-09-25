#!/usr/bin/env nextflow

include { UNTAR_CR } from './modules/local/tar/main.nf'
include { DOUBLETFINDER } from './modules/local/doubletFinder/main.nf'
include { SOUPX } from './modules/local/soupX/main.nf'
include { COLLATE_OUTPUTS } from './modules/local/collate_outputs/main.nf'
include { CREATE_INITIAL_SEURAT } from './modules/local/create_initial_seurat/main.nf'
include { ANALYZE_SEURAT_OBJECT } from './modules/local/analyze_seurat_object/main.nf'
include { CREATE_IMAGES_DGE } from './modules/local/create_images_dge/main.nf'
include { COLLATE_ANALYSIS } from './modules/local/collate_anaylsis/main.nf'

def validate_inputs(param_obj){
    // single value possibilities
    def required_options = [
        organism: ["mouse", "human"],
        starting_data:["cellranger", "matrix"],
        soupx_start: ["out", "no_clusters", "h5"],
        split_layers_by: ["Sample", "Experiment"],
        scale_data_features: ["all", "variable"]
    ]
    // multi value possibilities
    def required_multi_options = [
        normalization_config: ["sct", "standard"],
        integration_config: ["cca", "harmony", "rpca"],
        visualization: ["feature", "violin", "ridge", "dot"]
    ]
    param_obj.each { k, v ->
        if (required_options.containsKey(k)){
            if (!required_options[k].contains(v)){
                error("Invalid option for parameter ${k}: ${v}. Valid options are: ${required_options[k]}")
            }
        }
        else if (required_multi_options.containsKey(k)){
            def vals = v.split(",")
            vals.each { val ->
                if (!required_multi_options[k].contains(val)){
                    error("Invalid option for parameter ${k}: ${val}. Valid options are: ${required_multi_options[k]}")
                }
            }
        }
    }
    // Validate optional dependent params
    def dependent_params = [
        cc_method: ["standard", "alternative"],
        azimuth_ref_human: ["adiposeref", "bonemarrowref", "fetusref", "heartref", "humancortexref","kidneyref", "lungref", "pancreasref", "pbmcref", "tonsilref"],
        azimuth_ref_mouse: ["mousecortexref"]
    ]
    if (param_obj.cc_regression == "y"){
        if (!dependent_params.cc_method.contains(param_obj.cc_method)){
            error("Invalid option for parameter cc_method: ${param_obj.cc_method}. Valid options are: ${dependent_params.cc_method}")
        }
    }
    if (param_obj.run_azimuth == "y"){
        if (param_obj.organism == "human" && !dependent_params.azimuth_ref_human.contains(param_obj.azimuth_ref)){
            error("Invalid option for parameter azimuth_ref: ${param_obj.azimuth_ref}. Valid options are: ${dependent_params.azimuth_ref_human}")
        }
        else if (param_obj.organism == "mouse" && !dependent_params.azimuth_ref_mouse.contains(param_obj.azimuth_ref)){
            error("Invalid option for parameter azimuth_ref: ${param_obj.azimuth_ref}. Valid options are: ${dependent_params.azimuth_ref_mouse}")
        }
    }
}

workflow {
    main:
    validate_inputs(params)
    sample_list = Channel.fromList(params.sample_list)
    condition_list = Channel.fromList(params.condition_list).collect()
    input_dir_list = params.input_dir_list ? Channel.fromPath(params.input_dir_list.class == String ? params.input_dir_list.split(',') as List : params.input_dir_list) : Channel.empty()
    input_cr_tar_list = params.input_cr_tar_list ? Channel.fromPath(params.input_cr_tar_list.class == String ? params.input_cr_tar_list.split(',') as List : params.input_cr_tar_list) : Channel.empty()
    // Create meta dict of common inputs to reduce param passing and to mimic snakemake yaml
    meta = [
        PROJECT: params.project,
        ORGANISM: params.organism,
        RPATH: params.r_lib_path,
        RUN_SOUPX: String.valueOf(!params.disable_soupx),
        SOUPX_START: params.soupx_start,
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
        cc_regression: params.cc_regression,
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
    // input_list.view()
    UNTAR_CR(
        input_cr_tar_list
    )
    // sample IDs baked into dir structures from tar files
    // Will leverage that to create per-tar sample list

    paired_sample_dir = UNTAR_CR.out.flatMap { sample_str, dir_list ->
        if (dir_list instanceof List) {
            return [sample_str.tokenize("\n"), dir_list].transpose()
        } else {
            return [sample_str.tokenize("\n"), [dir_list]].transpose()
        }
    }

    paired_sample_dir.view()

    // if (!params.disable_doubletfinder){
    //     DOUBLETFINDER(
    //         meta,
    //         input_list
    //     )
    // }
    // if (!params.disable_soupx){
    //     SOUPX(
    //         meta,
    //         input_list
    //     )
    // }
    // // collate results so that next step can use them all together
    // samples = SOUPX.out.map { it[0] }.concat(DOUBLETFINDER.out.map { it[0] }).collect()
    // input_dirs = SOUPX.out.map { it[1] }.collect().combine(DOUBLETFINDER.out.map { it[1] }.collect())
    // COLLATE_OUTPUTS(
    //     meta,
    //     samples,
    //     input_dirs
    // )
    // seurat_filename = "data/endpoints/$params.project/analysis/RDS/${params.project}_initial_seurat_object.qs"
    // sample_list_flat = sample_list.collect()
    // input_dir_list_flat = input_dir_list.collect()
    // CREATE_INITIAL_SEURAT(
    //     COLLATE_OUTPUTS.out,
    //     sample_list_flat,
    //     condition_list,
    //     input_dir_list_flat,
    //     seurat_filename
    // )
    // ANALYZE_SEURAT_OBJECT(
    //     CREATE_INITIAL_SEURAT.out.seurat_file,
    //     params.aso_memory
    // )
    // CREATE_IMAGES_DGE(
    //     params.storage,
    //     ANALYZE_SEURAT_OBJECT.out.analyzed_seurat_object_file
    // )
    // analysis_dirs = CREATE_INITIAL_SEURAT.out.analysis_dir.combine(ANALYZE_SEURAT_OBJECT.out.analysis_path).combine(CREATE_IMAGES_DGE.out)
    // COLLATE_ANALYSIS(
    //     analysis_dirs
    // )
}