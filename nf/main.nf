#!/usr/bin/env nextflow

include { DOUBLETFINDER } from './modules/local/doubletFinder/main.nf'

workflow {
    main:
    sample = Channel.value(params.sample)
    starting_data = Channel.value(params.starting_data)
    input_dir = Channel.fromPath(params.input_dir)
    output_dir = map {}

    DOUBLETFINDER(
        sample,
        starting_data.collect()
        input_dir,
        output_dir
        params.mito_fraction.collect()
        params.min_feature_threshold.collect()
        params.int_components.collect()
        params.organism.collect()
    )
}