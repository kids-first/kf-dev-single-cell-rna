#!/usr/bin/env nextflow

include { DOUBLETFINDER } from '../../../modules/local/doubletFinder/main.nf'
include { SOUPX } from '../../../modules/local/soupX/main.nf'
include { COLLATE_OUTPUTS } from '../../../modules/local/collate_outputs/main.nf'
include { TAR_OUTPUTS as TAR_OUTPUTS_DBL } from '../../../modules/local/tar/main.nf'
include { TAR_OUTPUTS as TAR_OUTPUTS_SOUP } from '../../../modules/local/tar/main.nf'


workflow data_cleanup {
    take:
        meta
        src_sample_dir
    main:

    // if you've given matched doublet finder and soupX data, then you probably don't need to run doubletFinder again on that data
    if (!params.disable_doubletfinder){
        def doubletfinder_samples = src_sample_dir.doubletfinder.map { it[1] }.collect()
        dbl_input = src_sample_dir.cellranger
            .concat(src_sample_dir.matrix.filter {  !(doubletfinder_samples.contains(it[1])) })
            .map { src, sample, _condition, dir -> [src, sample, dir] }
        DOUBLETFINDER(
            meta,
            dbl_input
        )
        def dbl_tar_input = DOUBLETFINDER.out.map {_sample, dir -> ["", dir] }
        TAR_OUTPUTS_DBL(
            dbl_tar_input
        )
    }
    if (!params.disable_soupx){
        soupx_input = src_sample_dir.cellranger.map { src, sample, _condition, dir -> [src, sample, dir] }
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
        .concat(src_sample_dir.doubletfinder.map { _src, sample, _condition, _dir -> sample })
        .concat(src_sample_dir.matrix.map { _src, sample, _condition, _dir -> sample })
        .collect()
    input_dirs = SOUPX.out.map { _sample, dir -> dir }
        .concat(DOUBLETFINDER.out.map { _sample, dir -> dir })
        .concat(src_sample_dir.doubletfinder.map {  _src, _sample, _condition, dir -> dir })
        .concat(src_sample_dir.matrix.map { _src, _sample, _condition, dir -> dir })
        .collect()
    input_dirs.view()
    COLLATE_OUTPUTS(
        meta,
        samples,
        input_dirs
    )
}