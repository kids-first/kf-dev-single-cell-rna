#!/usr/bin/env nextflow

include { COUNT } from './modules/local/cell_ranger/count/main.nf'

workflow {
    main:
    sample = params.sample
    reads = params.reads ? Channel.fromPath(params.reads.class == String ? params.reads.split(',') as List : params.reads).collect() : Channel.empty()
    mates = params.mates ? Channel.fromPath(params.mates.class == String ? params.mates.split(',') as List : params.mates).collect() : Channel.empty()
    indices = params.indices ? Channel.fromPath(params.indices.class == String ? params.indices.split(',') as List : params.indices).collect() : Channel.value([])
    transcriptome = file(params.transcriptome)
    create_bam = String.valueOf(params.create_bam)

    indices ? indices.collect() : Channel.value([])
    COUNT(sample, create_bam, reads, mates, transcriptome, indices)
}