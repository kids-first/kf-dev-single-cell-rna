#!/usr/bin/env nextflow

include { DOUBLETFINDER } from './modules/local/doubletFinder/main.nf'
include { SOUPX } from './modules/local/soupX/main.nf'

workflow {
    main:
    sample = Channel.fromList(params.sample).view()
    project = params.project
    starting_data = Channel.value(params.starting_data)
    input_dir = Channel.fromPath(params.input_dir, type: 'dir').view()
    output_base = sample.map {"data/endpoints/$project/$it"}
    doubletFinder_output_dir = output_base.map { dir -> "$dir/doubletFinder" }
    soupx_output_dir = output_base.map { dir -> "$dir/soupX" }
    doubletFinder_output_dir.view()

    DOUBLETFINDER(
        sample,
        starting_data,
        input_dir,
        doubletFinder_output_dir,
        params.mito_cutoff,
        params.min_feature_threshold,
        params.int_components,
        params.organism
    )
    SOUPX(
        sample,
        params.data_type,
        input_dir,
        soupx_output_dir    )
}