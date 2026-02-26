#!/usr/bin/env nextflow

include { ANALYZE_SEURAT_OBJECT } from './modules/local/analyze_seurat_object/main.nf'
include { CREATE_IMAGES_DGE } from './modules/local/create_images_dge/main.nf'
include { COLLATE_ANALYSIS } from './modules/local/collate_anaylsis/main.nf'
include { TAR_OUTPUTS } from './modules/local/tar/main.nf'


def validate_inputs(param_obj){
    // single value possibilities
    def required_options = [
        organism: ["mouse", "human"],
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
        azimuth_ref_human: ["adiposeref", "bonemarrowref", "fetusref", "heartref", "humancortexref","kidneyref", "lungref", "pancreasref", "pbmcref", "tonsilref"],
        azimuth_ref_mouse: ["mousecortexref"]
    ]
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
    initial_seurat_object = channel.fromPath(params.initial_seurat_object) // TSV file with required headers: sample_id, condition, name, input_type. Optional header remap
    validate_inputs(params)
    meta = channel.value(
        [
            PROJECT: params.project,
            ORGANISM: params.organism,
            RPATH: params.r_lib_path,
            MITO: params.mito_cutoff,
            RIBO: params.ribo_cutoff,
            SEURAT_NORMALIZATION_METHOD: params.normalization_config,
            SEURAT_INTEGRATION_METHOD: params.integration_config,
            RESOLUTION: params.resolution_config,
            COMPONENTS: params.int_components,
            MITO_REGRESSION: params.mito_regression,
            RIBO_REGRESSION: params.ribo_regression,
            CELL_CYCLE_REGRESSION: params.cc_method ? 'y' : 'n',
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
    )

    analyze_seurat_object_input = meta.combine(initial_seurat_object)
    ANALYZE_SEURAT_OBJECT(
        analyze_seurat_object_input,
        params.aso_memory
    )

    CREATE_IMAGES_DGE(
        params.storage,
        ANALYZE_SEURAT_OBJECT.out.analyzed_seurat_object_file
    )
    analysis_dirs = ANALYZE_SEURAT_OBJECT.out.analysis_path.combine(CREATE_IMAGES_DGE.out)
    COLLATE_ANALYSIS(
        analysis_dirs
    )
    tar_output_input = COLLATE_ANALYSIS.out.map { dir_path -> ["data/endpoints/$params.project", dir_path] }
    TAR_OUTPUTS(
        tar_output_input
    )
}