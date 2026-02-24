#!/usr/bin/env nextflow

include { format_inputs } from './subworkflows/local/format_inputs/main.nf'
include { data_cleanup } from './subworkflows/local/data_cleanup/main.nf'
include { run_qc } from './subworkflows/local/run_qc/main.nf'

def validate_manifest(manifest){
    // single value possibilities
    def required_free_text = [
        "sample_id",
        "condition",
        "name",
    ]
    // multi value possibilities (lists)
    def required_enum = [
        input_type: [
            "h5_raw",
            "h5_filtered",
            "dir_cellranger",
            "tar_cellranger_count",
            "tar_cellranger_multi",
            "dir_doubletFinder",
            "tar_doubletFinder",
            "dir_soupX",
            "tar_soupX"
        ]
    ]
    manifest.splitCsv(header: true, sep: "\t").map {
        row -> row.each{
            required_free_text.each { required_key ->
                if (!row.containsKey(required_key) || row[required_key] == null || row[required_key].isEmpty()){
                    error("Manifest is missing required key ${required_key} or it is empty in row: ${row}")
                }
            }
            required_enum.each { required_key, valid_options ->
                if (!row.containsKey(required_key) || row[required_key] == null || row[required_key].isEmpty()){
                    error("Manifest is missing required key ${required_key} or it is empty in row: ${row}")
                }
                else {
                    def value = row[required_key]
                    if (!valid_options.contains(value)){
                        error("Invalid option for parameter ${required_key}: ${value}. Valid options are: ${valid_options}")
                    }
                }
            }
        }
    }

}

workflow {
    main:
    input_sample_sheet = channel.fromPath(params.input_sample_sheet) // TSV file with required headers: sample_id, condition, name, input_type. Optional header remap
    validate_manifest(input_sample_sheet)

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