#!/usr/bin/env nextflow

include { COUNT } from './modules/local/cell_ranger/count/main.nf'
include { MULTI } from './modules/local/cell_ranger/multi/main.nf'

workflow {
    main:
    // general
    mode = params.mode
    reads = params.reads ? Channel.fromPath(params.reads.class == String ? params.reads.split(',') as List : params.reads).collect() : Channel.empty()
    mates = params.mates ? Channel.fromPath(params.mates.class == String ? params.mates.split(',') as List : params.mates).collect() : Channel.empty()
    transcriptome = file(params.transcriptome)
    create_bam = String.valueOf(params.create_bam)
    // count specific
    sample = params.sample
    indices = params.indices ? Channel.fromPath(params.indices.class == String ? params.indices.split(',') as List : params.indices).collect() : Channel.value([])
    indices ? indices.collect() : Channel.value([])
    // multi specific
    sample_sheet = params.sample_csv ? Channel.fromPath(params.sample_csv) : Channel.empty()
    probe_set = params.probe_set ? Channel.fromPath(params.probe_set) : Channel.empty()
    library_fastq_id = params.library_fastq_id
    feature_types = params.feature_types

    if (mode == "count"){
        COUNT(
            sample,
            create_bam,
            reads,
            mates,
            transcriptome,
            indices
        )
    }

    if (mode == "multi"){
        MULTI(
            library_fastq_id,
            create_bam,
            reads,
            mates,
            transcriptome,
            feature_types,
            sample_sheet,
            probe_set
        )
    }
}