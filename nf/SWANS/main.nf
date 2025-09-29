#!/usr/bin/env nextflow

include { UNTAR_CR } from './modules/local/untar/main.nf'
include { TAR_OUTPUTS } from './modules/local/tar/main.nf'
include { TAR_OUTPUTS as TAR_OUTPUTS_DBL } from './modules/local/tar/main.nf'
include { TAR_OUTPUTS as TAR_OUTPUTS_SOUP } from './modules/local/tar/main.nf'
include { DOUBLETFINDER } from './modules/local/doubletFinder/main.nf'
include { SOUPX } from './modules/local/soupX/main.nf'
include { COLLATE_OUTPUTS } from './modules/local/collate_outputs/main.nf'
include { CREATE_INITIAL_SEURAT } from './modules/local/create_initial_seurat/main.nf'
include { ANALYZE_SEURAT_OBJECT } from './modules/local/analyze_seurat_object/main.nf'
include { CREATE_IMAGES_DGE } from './modules/local/create_images_dge/main.nf'
include { COLLATE_ANALYSIS } from './modules/local/collate_anaylsis/main.nf'


def parse_map_file(file_text){
    def sample_map = [:]
    def lines = file_text.split("\n")
    lines.drop(1).each { line ->
        def (sample, condition, remap) = line.tokenize("\t")
        if (remap == ""){
            remap = sample
        }
        sample_map[sample] = [condition, remap]
    }
    return sample_map
}
def validate_inputs(param_obj){
    // single value possibilities
    def required_options = [
        organism: ["mouse", "human"],
        soupx_start: ["out", "no_clusters", "h5"],
        split_layers_by: ["Sample", "Experiment"],
        scale_data_features: ["all", "variable"]
    ]
    // multi value possibilities
    def required_multi_options = [
        normalization_config: ["sct", "standard"],
        integration_config: ["cca", "harmony", "rpca"],
        visualization: ["feature", "violin", "ridge", "dot"]
    ]
    param_obj.each { k, v ->
        if (required_options.containsKey(k)){
            if (!required_options[k].contains(v)){
                error("Invalid option for parameter ${k}: ${v}. Valid options are: ${required_options[k]}")
            }
        }
        else if (required_multi_options.containsKey(k)){
            def vals = v.split(",")
            vals.each { val ->
                if (!required_multi_options[k].contains(val)){
                    error("Invalid option for parameter ${k}: ${val}. Valid options are: ${required_multi_options[k]}")
                }
            }
        }
    }
    // Validate optional dependent params
    def dependent_params = [
        cc_method: ["standard", "alternative"],
        azimuth_ref_human: ["adiposeref", "bonemarrowref", "fetusref", "heartref", "humancortexref","kidneyref", "lungref", "pancreasref", "pbmcref", "tonsilref"],
        azimuth_ref_mouse: ["mousecortexref"]
    ]
    if (param_obj.cc_regression == "y"){
        if (!dependent_params.cc_method.contains(param_obj.cc_method)){
            error("Invalid option for parameter cc_method: ${param_obj.cc_method}. Valid options are: ${dependent_params.cc_method}")
        }
    }
    if (param_obj.run_azimuth == "y"){
        if (param_obj.organism == "human" && !dependent_params.azimuth_ref_human.contains(param_obj.azimuth_ref)){
            error("Invalid option for parameter azimuth_ref: ${param_obj.azimuth_ref}. Valid options are: ${dependent_params.azimuth_ref_human}")
        }
        else if (param_obj.organism == "mouse" && !dependent_params.azimuth_ref_mouse.contains(param_obj.azimuth_ref)){
            error("Invalid option for parameter azimuth_ref: ${param_obj.azimuth_ref}. Valid options are: ${dependent_params.azimuth_ref_mouse}")
        }
    }
}

def parse_input_dir_src(dir_channel, src_channel, sample_map){
    // collate src and dir, parse out sample name from dir, and assign to desired name (could be same as original)
    return src_channel.merge(dir_channel).map { src, dir -> 
        if (src.toLowerCase() == "doubletfinder") {
            def sname_matcher = dir =~ /([^\/]+)[_\/]doubletFinder/
            def preceding_str = sname_matcher ? sname_matcher[0][1] : error("Could not parse sample name from input_dir_list entry ${dir}. Ensure directory name contains 'doubletFinder'.")
            return [src, sample_map[preceding_str][1], sample_map[preceding_str][0], dir]
        }
        else if (src.toLowerCase() == "soupx") {
            def sname_matcher = dir =~ /([^\/]+)[_\/]soupX/
            def preceding_str = sname_matcher ? sname_matcher[0][1] : error("Could not parse sample name from input_dir_list entry ${dir}. Ensure directory name contains 'soupX'.")
            return [src, sample_map[preceding_str][1], sample_map[preceding_str][0], dir]
        }
        else {
            return [src, sample_map[dir.name][1], sample_map[dir.name][0], dir]
        }
    }
}

def process_untar_outputs(untar_output, sample_map){
    return untar_output.flatMap { data_src, sample_str, dir_list ->
        def sample_list = sample_str.tokenize("\n")
        def remap = []
        def cond = []
        sample_list.each { sname ->
            if (sample_map.containsKey(sname)){
                remap << sample_map[sname][1] ?: sname
                cond << sample_map[sname][0]
            } else {
                error("Sample name ${sname} not found in sample_condition_map_file. Please ensure all samples are mapped.")
            }
        }
        if (dir_list instanceof List) {
            return [[data_src] * dir_list.size(), remap, cond, dir_list].transpose()
        } else {
            return [[data_src], remap, cond, [dir_list]].transpose()
        }
    }
}
workflow {
    main:
    validate_inputs(params)
    sample_condition_map_file = file(params.sample_condition_map_file)
    input_dir_list = params.input_dir_list ? Channel.fromPath(params.input_dir_list.class == String ? params.input_dir_list.split(',') as List : params.input_dir_list) : Channel.empty()
    input_dir_src_list = params.input_dir_src_list ? Channel.fromList(params.input_dir_src_list) : Channel.value([])
    input_tar_list = params.input_tar_list ? Channel.fromPath(params.input_tar_list.class == String ? params.input_tar_list.split(',') as List : params.input_tar_list) : Channel.empty()
    input_tar_src_list = params.input_tar_src_list ? Channel.fromList(params.input_tar_src_list) : Channel.value([])
    // Create meta dict of common inputs to reduce param passing and to mimic snakemake yaml
    meta = [
        PROJECT: params.project,
        ORGANISM: params.organism,
        RPATH: params.r_lib_path,
        RUN_SOUPX: String.valueOf(!params.disable_soupx),
        SOUPX_START: params.soupx_start,
        RUN_DOUBLETFINDER: String.valueOf(!params.disable_doubletfinder),
        MITO: params.mito_cutoff,
        RIBO: params.ribo_cutoff,
        MIN_FEATURE_THRESHOLD: params.min_feature_threshold,
        MAX_FEATURE_THRESHOLD: params.max_feature_threshold,
        SEURAT_NORMALIZATION_METHOD: params.normalization_config,
        SEURAT_INTEGRATION_METHOD: params.integration_config,
        RESOLUTION: params.resolution_config,
        COMPONENTS: params.int_components,
        MITO_REGRESSION: params.mito_regression,
        RIBO_REGRESSION: params.ribo_regression,
        cc_regression: params.cc_regression,
        NUM_VARIABLE_FEATURES: params.num_var_features,
        SCALE_DATA_FEATURES: params.scale_data_features,
        SPLIT_LAYERS_BY: params.split_layers_by,
        REFERENCE_BASED_INTEGRATION: params.ref_based_integration,
        RUN_AZIMUTH: params.run_azimuth,
        RUN_TRANSFERDATA: params.run_transferdata,
        TSNE: params.include_tsne,
        CONSERVED_GENES: params.conserved_genes,
        VISUALIZATION: params.visualization
    ]
    // dir names typically drive sample names, but not always desired. use sample map to enforce desired names
    input_meta_tar = input_tar_src_list.merge(input_tar_list).map { src, tar -> [src, tar] }
    UNTAR_CR(
        input_meta_tar
    )
    // parse sample map file
    sample_condition_map =  parse_map_file(sample_condition_map_file.text)
    print(sample_condition_map)
    // initialize dir, tar, and src data as one channel with src, sample_name, condition, dir
    parse_input_dir_src(input_dir_list, input_dir_src_list, sample_condition_map).concat(process_untar_outputs(UNTAR_CR.out, sample_condition_map)).branch{ parsed_input -> 
        doubletfinder: parsed_input[0].toLowerCase() == "doubletfinder"
        matrix: parsed_input[0].toLowerCase() == "soupx"
        cellranger: parsed_input[0].toLowerCase() == "cellranger"
    }.set{src_sample_dir}


    if (!params.disable_doubletfinder){
        dbl_input = src_sample_dir.cellranger.concat(src_sample_dir.matrix).map { src, sample, _condition, dir -> [src, sample, dir] }
        DOUBLETFINDER(
            meta,
            dbl_input
        )
        def dbl_tar_input = DOUBLETFINDER.out.map {_sample, dir -> ["", dir] }
        TAR_OUTPUTS_DBL(
            dbl_tar_input
        )
    }
    if (!params.disable_soupx){
        soupx_input = src_sample_dir.cellranger.map { src, sample, _condition, dir -> [src, sample, dir] }
        SOUPX(
            meta,
            soupx_input
        )
        def soupx_tar_input = SOUPX.out.map { _sample, dir -> ["", dir] }
        TAR_OUTPUTS_SOUP(
            soupx_tar_input
        )
    }
    // // collate results so that next step can use them all together
    samples = SOUPX.out.map { it[0] }.concat(DOUBLETFINDER.out.map { it[0] }).collect()
    input_dirs = SOUPX.out.map { it[1] }.collect().combine(DOUBLETFINDER.out.map { it[1] }.collect())
    COLLATE_OUTPUTS(
        meta,
        samples,
        input_dirs
    )
    seurat_filename = "data/endpoints/$params.project/analysis/RDS/${params.project}_initial_seurat_object.qs"
    // use metadata from matrix and cellranger from src_sample_dir to help collate create the initial sample list file
    def (sample_list_flat, condition_list, input_dir_list_flat) = [ [],  [], [] ]
    src_sample_dir.cellranger.concat(src_sample_dir.matrix).map { _src, sample, condition, dir ->
        sample_list_flat << sample
        condition_list << condition
        input_dir_list_flat << dir
    }
    CREATE_INITIAL_SEURAT(
        COLLATE_OUTPUTS.out,
        sample_list_flat,
        condition_list,
        input_dir_list_flat,
        seurat_filename
    )
    ANALYZE_SEURAT_OBJECT(
        CREATE_INITIAL_SEURAT.out.seurat_file,
        params.aso_memory
    )
    CREATE_IMAGES_DGE(
        params.storage,
        ANALYZE_SEURAT_OBJECT.out.analyzed_seurat_object_file
    )
    analysis_dirs = CREATE_INITIAL_SEURAT.out.analysis_dir.combine(ANALYZE_SEURAT_OBJECT.out.analysis_path).combine(CREATE_IMAGES_DGE.out)
    COLLATE_ANALYSIS(
        analysis_dirs
    )
    TAR_OUTPUTS(
        COLLATE_ANALYSIS.out.map { ["data/endpoints/$params.project", it[1]] }
    )
}