#!/usr/bin/env nextflow

include { DOUBLETFINDER } from './modules/local/doubletFinder/main.nf'
include { SOUPX } from './modules/local/soupX/main.nf'

workflow {
    main:
    sample_list = Channel.fromList(params.sample_list)
    project = params.project
    starting_data = Channel.value(params.starting_data)
    input_dir_list = Channel.fromPath(params.input_dir_list, type: 'dir')
    output_base = sample_list.map {"data/endpoints/$project/$it"}
    doubletfinder_output_dir = output_base.map { dir -> "$dir/doubletFinder" }
    soupx_output_dir = output_base.map { dir -> "$dir/soupX/" }
    doubletfinder_output_dir

    input_list = sample_list.merge(input_dir_list).map { sample, input_dir ->
        [sample, input_dir]
    }.view()
    if (!params.disable_soupx){
        SOUPX(
            input_list,
            params.data_type,
            soupx_output_dir,
            starting_data
        )
        SOUPX.out.filtered_counts_dir
    }
    // if (!params.disable_doubletfinder){
    //     DOUBLETFINDER(
    //         input_list,
    //         starting_data,
    //         doubletfinder_output_dir,
    //         params.mito_cutoff,
    //         params.min_feature_threshold,
    //         params.int_components,
    //         params.organism
    //     )
    // }

}