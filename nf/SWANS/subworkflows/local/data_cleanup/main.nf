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
        doubletfinder_samples = doublet_data_dir.map { metadata, _dir -> metadata.sample }.collect()
        dbl_input = cellranger_data_dir
            .concat(matrix_data_dir.filter { metadata, _dir -> !(doubletfinder_samples.contains(metadata.sample)) })
        DOUBLETFINDER(
            meta,
            dbl_input
        )
        dbl_tar_input = DOUBLETFINDER.out.map {_metadata, dir -> ["", dir] }
        TAR_OUTPUTS_DBL(
            dbl_tar_input
        )
    }
    if (!params.disable_soupx){
        SOUPX(
            meta,
            cellranger_data_dir
        )
        def soupx_tar_input = SOUPX.out.map { _metadata, dir -> ["", dir] }
        samples = SOUPX.out.map { meta_data, _dir -> meta_data.sample_id }
        input_dirs = SOUPX.out.map { _metadata, dir -> dir }

        TAR_OUTPUTS_SOUP(
            soupx_tar_input
        )
    }
    // COLLATE RESULTS
    // initialize with cellranger results if soupX disabled
    if (params.disable_soupx){
        samples = cellranger_data_dir.map { meta_data, _dir -> meta_data.sample_id }
        input_dirs = cellranger_data_dir.map { _metadata, dir -> dir }
    }
    if (!params.disable_doubletfinder) {
        samples = samples.concat(DOUBLETFINDER.out.map { meta_data, _dir -> meta_data.sample_id })
        input_dirs = input_dirs.concat(DOUBLETFINDER.out.map { _metadata, dir -> dir })
    }
    samples = samples.concat(doublet_data_dir.out.map { meta_data, _dir -> meta_data.sample_id })
        .concat(matrix_data_dir.out.map { meta_data, _dir -> meta_data.sample_id })
        .collect()
    input_dirs = input_dirs.concat(doublet_data_dir.map { _metadata, dir -> dir })
        .concat(matrix_data_dir.map { _metadata, dir -> dir })
        .collect()
    COLLATE_OUTPUTS(
        meta,
        samples,
        input_dirs
    )

    emit:
    COLLATE_OUTPUTS.out
}