#!/usr/bin/env nextflow

include { COUNT } from './modules/local/cell_ranger/count/main.nf'
include { MULTI } from './modules/local/cell_ranger/multi/main.nf'
include { UNTAR_REF } from './modules/local/untar_ref/main.nf'
include { TAR_DIR } from './modules/local/tar_dir/main.nf'

workflow {
    main:
    // general
    reads = params.reads ? Channel.fromPath(params.reads.class == String ? params.reads.split(',') as List : params.reads).collect() : Channel.empty()
    mates = params.mates ? Channel.fromPath(params.mates.class == String ? params.mates.split(',') as List : params.mates).collect() : Channel.empty()
    transcriptome_dir = params.transcriptome_dir ? Channel.fromPath(params.transcriptome_dir) : ""
    transcriptome_tar = params.transcriptome_tar ? Channel.fromPath(params.transcriptome_tar) : ""
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
            reads,
            mates,
            transcriptome_dir,
            indices
        )
        sample_analysis_dir = sample.merge(COUNT.out.analysis_dir)
        TAR_DIR(sample_analysis_dir)
    }

    else{
        MULTI(
            library_fastq_id,
            reads,
            mates,
            transcriptome_dir,
            feature_types,
            sample_sheet,
            probe_set
        )
    }
    pattern: /.*\/per_sample_outs/([^\/]+?)\/count\/analysis/
    multi_sample_analysis_dir = MULTI.out.multi_analysis.map{
        dirname ->

    }

}