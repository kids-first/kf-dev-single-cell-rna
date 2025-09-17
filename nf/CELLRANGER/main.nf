#!/usr/bin/env nextflow

include { COUNT } from './modules/local/cell_ranger/count/main.nf'
include { MULTI } from './modules/local/cell_ranger/multi/main.nf'
include { UNTAR_REF } from './modules/local/tar/main.nf'

workflow {
    main:
    // general
    mode = params.mode
    reads = params.reads ? Channel.fromPath(params.reads.class == String ? params.reads.split(',') as List : params.reads).collect() : Channel.empty()
    mates = params.mates ? Channel.fromPath(params.mates.class == String ? params.mates.split(',') as List : params.mates).collect() : Channel.empty()
    transcriptome_dir = params.transcriptome_dir ? file(params.transcriptome_dir) : ""
    transcriptome_tar = params.transcriptome_tar ? file(params.transcriptome_tar) : ""
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

    if (transcriptome_dir == "" && transcriptome_tar != ""){
        transcriptome_dir = UNTAR_REF(transcriptome_tar)
    } else if ((transcriptome_dir == "" && transcriptome_tar == "") || (transcriptome_dir != "" && transcriptome_tar != "")){
        error "Must provide one of either a path to a transcriptome directory or a tar file of the reference!"
    }
    if (mode == "count"){
        COUNT(
            sample,
            create_bam,
            reads,
            mates,
            transcriptome_dir,
            indices
        )
    }

    if (mode == "multi"){
        MULTI(
            library_fastq_id,
            create_bam,
            reads,
            mates,
            transcriptome_dir,
            feature_types,
            sample_sheet,
            probe_set
        )
    }
}