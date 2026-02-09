#!/usr/bin/env nextflow

include { DOUBLETFINDER } from '../../../modules/local/doubletFinder/main.nf'
include { SOUPX } from '../../../modules/local/soupX/main.nf'
include { COLLATE_OUTPUTS } from '../../../modules/local/collate_outputs/main.nf'
include { TAR_OUTPUTS as TAR_OUTPUTS_DBL } from '../../../modules/local/tar/main.nf'
include { TAR_OUTPUTS as TAR_OUTPUTS_SOUP } from '../../../modules/local/tar/main.nf'


workflow data_cleanup {
    take:
        doublet_data_dir
        matrix_data_dir
        cellranger_data_dir
    main:
    // Create meta dict of common inputs to reduce param passing and to mimic snakemake yaml
    meta = [
        PROJECT: params.project,
        ORGANISM: params.organism,
        RPATH: params.r_lib_path,
        RUN_SOUPX: String.valueOf(!params.disable_soupx),
        SOUPX_START: params.soupx_start,
        RUN_DOUBLETFINDER: String.valueOf(!params.disable_doubletfinder),
        MITO: params.mito_cutoff,
        RIBO: params.ribo_cutoff,
        MIN_FEATURE_THRESHOLD: params.min_feature_threshold,
        MAX_FEATURE_THRESHOLD: params.max_feature_threshold,
        COMPONENTS: params.int_components
    ]

    // if you've given matched doublet finder and soupX data, then you probably don't need to run doubletFinder again on that data
    if (!params.disable_doubletfinder){
        doubletfinder_samples = doublet_data_dir.map { _src, sample, _condition, _dir -> sample }.collect()
        dbl_input = cellranger_data_dir
            .concat(matrix_data_dir.filter { _src, sample, _condition, _dir -> !(doubletfinder_samples.contains(sample)) })
            .map { src, sample, _condition, dir -> [src, sample, dir] }
        DOUBLETFINDER(
            meta,
            dbl_input
        )
        dbl_tar_input = DOUBLETFINDER.out.map {_sample, dir -> ["", dir] }
        TAR_OUTPUTS_DBL(
            dbl_tar_input
        )
    }
    if (!params.disable_soupx){
        soupx_input = cellranger_data_dir.map { src, sample, _condition, dir -> [src, sample, dir] }
        SOUPX(
            meta,
            soupx_input
        )
        def soupx_tar_input = SOUPX.out.map { _sample, dir -> ["", dir] }
        TAR_OUTPUTS_SOUP(
            soupx_tar_input
        )
    }
    // COLLATE RESULTS
    samples = SOUPX.out.map { sample, _dir -> sample }
        .concat(DOUBLETFINDER.out.map { sample, _dir -> sample })
        .concat(doublet_data_dir.map { _src, sample, _condition, _dir -> sample })
        .concat(matrix_data_dir.map { _src, sample, _condition, _dir -> sample })
        .collect()
    input_dirs = SOUPX.out.map { _sample, dir -> dir }
        .concat(DOUBLETFINDER.out.map { _sample, dir -> dir })
        .concat(doublet_data_dir.map {  _src, _sample, _condition, dir -> dir })
        .concat(matrix_data_dir.map { _src, _sample, _condition, dir -> dir })
        .collect()
    COLLATE_OUTPUTS(
        meta,
        samples,
        input_dirs
    )

    emit:
    COLLATE_OUTPUTS.out
}