#!/usr/bin/env nextflow

include { TAR_OUTPUTS } from '../../../modules/local/tar/main.nf'
include { CREATE_INITIAL_SEURAT } from '../../../modules/local/create_initial_seurat/main.nf'

workflow run_qc {
    take:
        cellranger_data_dir
        matrix_data_dir
        cleanup_dir
    main:
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
    tar_output_input = CREATE_INITIAL_SEURAT.out.analysis_dir.map { dir -> ["", dir] }
    TAR_OUTPUTS(
        tar_output_input
    )
    emit:
        seurat_qs = CREATE_INITIAL_SEURAT.out.seurat_qs
        qc_report = CREATE_INITIAL_SEURAT.out.qc_report
        seurat_tar = TAR_OUTPUTS.out

}