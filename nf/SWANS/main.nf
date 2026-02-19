#!/usr/bin/env nextflow

include { format_inputs } from './subworkflows/local/format_inputs/main.nf'
include { data_cleanup } from './subworkflows/local/data_cleanup/main.nf'
include { run_qc } from './subworkflows/local/run_qc/main.nf'

def validate_inputs(param_obj){
    // single value possibilities
    def valid_options = [
        organism: ["mouse", "human"],
        soupx_start: ["outs", "no_clusters", "h5"]
    ]
    // multi value possibilities (lists)
    def valid_multi_options = [
        input_dir_src_list: ["cellranger", "doubletFinder", "soupX"],
        input_file_src_list: ["h5_raw", "h5_filtered", "cellranger", "doubletFinder", "soupX"]
    ]

    param_obj.each { k, v ->
        if (valid_options.containsKey(k) && !(v == null || v.isEmpty())){
            if (!valid_options[k].contains(v)){
                error("Invalid option for parameter ${k}: ${v}. Valid options are: ${valid_options[k]}")
            }
        }
        else if (valid_multi_options.containsKey(k) && !(v == null || v.isEmpty())){
            def vals = v instanceof String ? v.split(",") : v
            vals.each { val ->
                if (!valid_multi_options[k].contains(val)){
                    error("Invalid option for parameter ${k}: ${val}. Valid options are: ${valid_multi_options[k]}")
                }
            }
        }
    }
}

workflow {
    main:
    validate_inputs(params)
    input_sample_sheet = channel.fromPath(params.input_sample_sheet) // TSV file with required headers: sample_id, condition, name, input_type. Optional header remap

    // FORMAT INPUTS
    (cellranger_data_dir, doublet_data_dir, matrix_data_dir) = format_inputs(
        input_sample_sheet
    )
    // weird quirk where if other matrices empty, values from one channel spill over?! need to initialize empty channels to avoid this
    cellranger_data_dir = (cellranger_data_dir ?: channel.value([]))
    doublet_data_dir = doublet_data_dir ?: channel.value([])
    matrix_data_dir = (matrix_data_dir ?: channel.value([]))
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