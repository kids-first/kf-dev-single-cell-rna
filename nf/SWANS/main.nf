#!/usr/bin/env nextflow

include { format_inputs } from './subworkflows/local/format_inputs/main.nf'
include { data_cleanup } from './subworkflows/local/data_cleanup/main.nf'
include { TAR_OUTPUTS } from './modules/local/tar/main.nf'
include { CREATE_INITIAL_SEURAT } from './modules/local/create_initial_seurat/main.nf'
include { ANALYZE_SEURAT_OBJECT } from './modules/local/analyze_seurat_object/main.nf'
include { CREATE_IMAGES_DGE } from './modules/local/create_images_dge/main.nf'
include { COLLATE_ANALYSIS } from './modules/local/collate_anaylsis/main.nf'


def validate_inputs(param_obj){
    // single value possibilities
    def required_options = [
        organism: ["mouse", "human"],
        soupx_start: ["outs", "no_clusters", "h5"],
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
    validate_inputs(params)
    sample_condition_map_file = file(params.sample_condition_map_file)
    input_dir_list = params.input_dir_list ? channel.fromPath(params.input_dir_list.class == String ? params.input_dir_list.split(',') as List : params.input_dir_list) : channel.empty()
    input_dir_src_list = params.input_dir_src_list ? channel.fromList(params.input_dir_src_list) : channel.value([])
    input_tar_list = params.input_tar_list ? channel.fromPath(params.input_tar_list.class == String ? params.input_tar_list.split(',') as List : params.input_tar_list) : channel.empty()
    input_tar_src_list = params.input_tar_src_list ? channel.fromList(params.input_tar_src_list) : channel.value([])

    // FORMAT INPUTS
    (meta, doublet_data_dir, matrix_data_dir, cellranger_data_dir) = format_inputs(
        input_tar_src_list,
        input_tar_list,
        sample_condition_map_file,
        input_dir_list,
        input_dir_src_list
    )

    // CLEAN UP DATA
    cleanup_dir = data_cleanup(
        meta,
        doublet_data_dir,
        matrix_data_dir,
        cellranger_data_dir
    )

    // CREATE SEURAT OBJ/QC
    seurat_filename = "data/endpoints/$params.project/analysis/RDS/${params.project}_initial_seurat_object.qs"
    // use metadata from matrix and cellranger from src_sample_dir to help collate create the initial sample list file
    def (sample_list_flat, condition_list, input_dir_list_flat) = [ [],  [], [] ]
    cellranger_data_dir.concat(matrix_data_dir).map { _src, sample, condition, dir ->
        sample_list_flat << sample
        condition_list << condition
        input_dir_list_flat << dir
    }
    CREATE_INITIAL_SEURAT(
        cleanup_dir,
        sample_list_flat,
        condition_list,
        input_dir_list_flat,
        seurat_filename
    )
    CREATE_INITIAL_SEURAT.out.seurat_file.view()
    CREATE_INITIAL_SEURAT.out.qc_report.view()
    tar_output_input = CREATE_INITIAL_SEURAT.out.analysis_dir.map { dir -> ["", dir] }
    TAR_OUTPUTS(
            tar_output_input
        )
}