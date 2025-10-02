#!/usr/bin/env nextflow

include { COUNT } from './modules/local/cell_ranger/count/main.nf'
include { MULTI } from './modules/local/cell_ranger/multi/main.nf'
include { UNTAR_REF } from './modules/local/tar/main.nf'

workflow {
    main:
    // general
    mode = Channel.value(params.mode)
    reads = params.reads ? Channel.fromPath(params.reads.class == String ? params.reads.split(',') as List : params.reads).collect() : Channel.empty()
    mates = params.mates ? Channel.fromPath(params.mates.class == String ? params.mates.split(',') as List : params.mates).collect() : Channel.empty()
    transcriptome_dir = params.transcriptome_dir ? Channel.fromPath(params.transcriptome_dir) : ""
    transcriptome_tar = params.transcriptome_tar ? Channel.fromPath(params.transcriptome_tar) : ""
    create_bam = String.valueOf(params.create_bam)
    // count specific
    sample = Channel.value(params.sample)
    indices = params.indices ? Channel.fromPath(params.indices.class == String ? params.indices.split(',') as List : params.indices).collect() : Channel.value([])
    // multi specific
    sample_sheet = params.sample_csv ? Channel.fromPath(params.sample_csv) : Channel.empty()
    probe_set = params.probe_set ? Channel.fromPath(params.probe_set) : Channel.empty()
    library_fastq_id = Channel.value(params.library_fastq_id)
    feature_types = Channel.value(params.feature_types)

    if (!params.transcriptome_tar && !params.transcriptome_dir) {
        error "Must provide one of either a path to a transcriptome directory or a tar file of the reference!"
    } else if (params.transcriptome_dir) {
        transcriptome_dir = Channel.fromPath(params.transcriptome_dir)
    } else {
        UNTAR_REF(transcriptome_tar)
        transcriptome_dir = UNTAR_REF.out
    }

    if (params.mode == "count"){
        COUNT(
            sample,
            create_bam,
            reads,
            mates,
            transcriptome_dir,
            indices
        )
    }

    else{
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