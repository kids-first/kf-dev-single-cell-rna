#!/usr/bin/env nextflow

include { format_inputs } from './subworkflows/local/format_inputs/main.nf'
include { data_cleanup } from './subworkflows/local/data_cleanup/main.nf'
include { run_qc } from './subworkflows/local/run_qc/main.nf'

def validate_inputs(param_obj){
    // single value possibilities
    def valid_options = [
        organism: ["mouse", "human"],
        soupx_start: ["outs", "no_clusters", "h5"],
        input_file_src_list: ["doubletFinder", "matrix", "cellranger", "h5_raw", "h5_filtered"],
        input_dir_src_list: ["doubletFinder", "matrix", "cellranger"]
    ]

    param_obj.each { k, v ->
        if (valid_options.containsKey(k) && !(v == null || v.isEmpty())){
            if (!valid_options[k].contains(v)){
                error("Invalid option for parameter ${k}: ${v}. Valid options are: ${valid_options[k]}")
            }
        }

    }
}

workflow {
    main:
    validate_inputs(params)
    sample_condition_map_file = file(params.sample_condition_map_file) // TSV file with header and three columns: sample, condition, remap. remap is optional and if not provided, sample name will be used as-is.
    input_dir_list = params.input_dir_list ? channel.fromPath(params.input_dir_list.class == String ? params.input_dir_list.split(',') as List : params.input_dir_list) : channel.empty() // input dirs from any of cell ranger or previous runs of doubletFinder or soupX. Required if no tar inputs
    input_dir_src_list = params.input_dir_src_list ? channel.fromList(params.input_dir_src_list) : channel.value([]) // list of sources corresponding to input_dir_list, like cellranger, doubletFinder, soupX. Required if no tar inputs. 
    input_file_list = params.input_file_list ? channel.fromPath(params.input_file_list.class == String ? params.input_file_list.split(',') as List : params.input_file_list) : channel.empty() // input files from cell ranger or previous runs of doubletFinder or soupX as tar balls or h5 raw + filtered. Required if no dir inputs
    input_file_src_list = params.input_file_src_list ? channel.fromList(params.input_file_src_list) : channel.value([]) // list of sources corresponding to input_file_list, like cellranger, doubletFinder, soupX. h5_raw, h5_filtered. Required if no dir inputs.

    // FORMAT INPUTS
    (doublet_data_dir, matrix_data_dir, cellranger_data_dir) = format_inputs(
        input_file_src_list,
        input_file_list,
        sample_condition_map_file,
        input_dir_list,
        input_dir_src_list
    )

    // CLEAN UP DATA
    cleanup_dir = data_cleanup(
        doublet_data_dir,
        matrix_data_dir,
        cellranger_data_dir
    )

    // CREATE SEURAT OBJ/QC
    run_qc(
        cellranger_data_dir,
        matrix_data_dir,
        cleanup_dir
    )
}