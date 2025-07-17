#!/usr/bin/env nextflow

include { DOUBLETFINDER } from './modules/local/doubletFinder/main.nf'

workflow {
    main:
    sample = Channel.fromList(params.sample).view()
    project = params.project
    starting_data = Channel.value(params.starting_data)
    input_dir = Channel.fromPath(params.input_dir, type: 'dir').view()
    output_dir = sample.map {"data/endpoints/$project/$it/doubletFinder"}
    output_dir.view()

    DOUBLETFINDER(
        sample,
        starting_data,
        input_dir,
        output_dir,
        params.mito_cutoff,
        params.min_feature_threshold,
        params.int_components,
        params.organism
    )
}