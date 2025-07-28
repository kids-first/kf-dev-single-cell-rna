#!/usr/bin/env nextflow

include { DOUBLETFINDER } from './modules/local/doubletFinder/main.nf'
include { SOUPX } from './modules/local/soupX/main.nf'
include { TEST } from './modules/local/test/main.nf'

workflow {
    main:
    sample_list = Channel.fromList(params.sample_list)
    starting_data = Channel.value(params.starting_data)
    input_dir_list = Channel.fromPath(params.input_dir_list, type: 'dir')
    data_dir = "data/endpoints"

    // input_list = sample_list.merge(input_dir_list).map { sample, input_dir ->
    //    [sample, input_dir]
    //}
    // if (!params.disable_doubletfinder){
    //     DOUBLETFINDER(
    //         sample_list,
    //         input_dir_list,
    //         starting_data,
    //         data_dir,
    //         params.mito_cutoff,
    //         params.min_feature_threshold,
    //         params.int_components,
    //         params.organism
    //     )
    // }
    if (!params.disable_soupx){
        SOUPX(
            sample_list,
            input_dir_list,
            params.data_type,
            data_dir,
            starting_data
        )
    }
    SOUPX.out.view()
    TEST(SOUPX.out)
    TEST.out.view()

}