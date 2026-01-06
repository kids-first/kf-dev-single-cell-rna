#!/usr/bin/env nextflow

include { COUNT } from './modules/local/cell_ranger/count/main.nf'
include { MULTI } from './modules/local/cell_ranger/multi/main.nf'
include { UNTAR_REF } from './modules/local/untar_ref/main.nf'
include { TAR_DIR } from './modules/local/tar_dir/main.nf'

def validate_inputs(param_obj){
    if (!param_obj.transcriptome_tar && !param_obj.transcriptome_dir) {
        error "Must provide one of either a path to a transcriptome directory or a tar file of the reference!"
    }
    def required_params = ["mode", "reads", "mates", "create_bam"]
    def required_count_params = ["sample"]
    def required_multi_params = ["library_fastq_id", "sample_csv", "feature_types"]
    def missing_params = []
    required_params.each{ param -> if (!param_obj[param]) missing_params << param }
    if (param_obj.mode == "count"){
        required_count_params.each{ param -> if (!param_obj[param]) missing_params << param }
    }
    if (param_obj.mode == "multi"){
        required_multi_params.each{ param -> if (!param_obj[param]) missing_params << param }
    }
    if (missing_params.size() > 0) {
        error "Missing required parameters for mode ${param_obj.mode}: ${missing_params.join(', ')}"
    }
}

workflow {
    main:
    validate_inputs(params)
    // general
    reads = params.reads ? channel.fromPath(params.reads.class == String ? params.reads.split(',') as List : params.reads).collect() : channel.empty()
    mates = params.mates ? channel.fromPath(params.mates.class == String ? params.mates.split(',') as List : params.mates).collect() : channel.empty()
    transcriptome_dir = params.transcriptome_dir ? channel.fromPath(params.transcriptome_dir) : ""
    transcriptome_tar = params.transcriptome_tar ? channel.fromPath(params.transcriptome_tar) : ""
    // count specific
    sample = channel.value(params.sample)
    indices = params.indices ? channel.fromPath(params.indices.class == String ? params.indices.split(',') as List : params.indices).collect() : channel.value([])
    // multi specific
    sample_sheet = params.sample_csv ? channel.fromPath(params.sample_csv) : channel.empty()
    probe_set = params.probe_set ? channel.fromPath(params.probe_set) : channel.empty()
    library_fastq_id = channel.value(params.library_fastq_id)
    feature_types = channel.value(params.feature_types)

    if (params.transcriptome_dir) {
        transcriptome_dir = channel.fromPath(params.transcriptome_dir)
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
        sample_analysis_dir = COUNT.out.analysis_dir.map { dirname -> [params.sample + ".cellranger_analysis", dirname] }
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
        // extract sample_id from path
        pattern = /.*\/per_sample_outs\/(\S+)\/count\/analysis/
        multi_sample_analysis_dir = MULTI.out.multi_analysis.flatten().map{
            dirname ->
            def sample_id_match = dirname =~ pattern
            return sample_id_match ? [sample_id_match[0][1] + ".cellranger_analysis", dirname] : error("Could not find sample id in path: ${dirname}")
        }
        TAR_DIR(multi_sample_analysis_dir)
    }

}