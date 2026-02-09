#!/usr/bin/env nextflow

include { UNTAR_CR } from '../../../modules/local/untar/main.nf'

def parse_map_file(file_text){
    def sample_map = [:]
    def lines = file_text.split("\n")
    lines.drop(1).each { line ->
        def (sample, condition, remap) = line.tokenize("\t")
        sample_map[sample] = [condition, remap ?: sample]
    }
    return sample_map
}

def process_untar_outputs(untar_output, sample_map, pattern_map){
    return untar_output.flatMap { data_src, sample_str, dir_list ->
        def sample_list = sample_str.tokenize("\n")
        def remap = []
        def cond = []
        sample_list.each { sname ->
            // extract sample name if coming from doubletfinder or soupx
            if (pattern_map[data_src.toLowerCase()]){
                def sname_matcher = sname =~ pattern_map[data_src.toLowerCase()]
                sname = sname_matcher ? sname_matcher[0][1] : error("Could not parse sample name from untar_outputs entry $sname. Ensure directory name contains $data_src.")
            }
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

def parse_input_dir_src(dir_channel, src_channel, sample_map, pattern_map){
    // collate src and dir, parse out sample name from dir, and assign to desired name (could be same as original)
    return src_channel.merge(dir_channel).map { src, dir -> 
        if (pattern_map[src.toLowerCase()]) {
            def sname_matcher = dir =~ pattern_map[src.toLowerCase()]
            def preceding_str = sname_matcher ? sname_matcher[0][1] : error("Could not parse sample name from input_dir_list entry $dir. Ensure directory name contains $src.")
            return [src, sample_map[preceding_str][1], sample_map[preceding_str][0], dir]
        }
        else {
            return [src, sample_map[dir.name][1], sample_map[dir.name][0], dir]
        }
    }
}

workflow format_inputs {
    take:
        input_tar_src_list
        input_tar_list
        sample_condition_map_file
        input_dir_list
        input_dir_src_list
    main:
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
        COMPONENTS: params.int_components
    ]
    // dir names typically drive sample names, but not always desired. use sample map to enforce desired names
    input_meta_tar = input_tar_src_list.merge(input_tar_list) { src, tar -> [src, tar] }
    UNTAR_CR(
        input_meta_tar
    )
    // parse sample map file
    sample_condition_map =  parse_map_file(sample_condition_map_file.text)
    print(sample_condition_map)
    // initialize dir, tar, and src data as one channel with src, sample_name, condition, dir
    dirname_pattern = [
        doubletfinder: /([^\/]+)[_\/]doubletFinder/,
        matrix: /([^\/]+)[_\/]soupX/
    ]
    parse_input_dir_src(input_dir_list, input_dir_src_list, sample_condition_map, dirname_pattern).concat(process_untar_outputs(UNTAR_CR.out, sample_condition_map, dirname_pattern)).branch{ parsed_input -> 
            doubletfinder: parsed_input[0].toLowerCase() == "doubletfinder"
            matrix: parsed_input[0].toLowerCase() == "matrix"
            cellranger: parsed_input[0].toLowerCase() == "cellranger"
        }.set{src_sample_dir}
    emit:
        meta
        doubletfinder = src_sample_dir.doubletfinder
        matrix = src_sample_dir.matrix 
        cellranger = src_sample_dir.cellranger
}