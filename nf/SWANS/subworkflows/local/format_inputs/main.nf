#!/usr/bin/env nextflow

include { UNTAR_CR } from '../../../modules/local/untar/main.nf'
include { CONVERT_H5 } from '../../../modules/local/convert_h5/main.nf'

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
    // takes output of UNTAR_CR, sample map, and pattern map, and parses to create a unified channel of src, sample, condition, dir.
    // pattern map used to parse sample name from dir name if coming from doubletFinder or soupX /(<sample_name>_<doubletfinder/soupx> in the dir name), otherwise assumed to be cell ranger output and parsed as-is
    // sample map used to assign desired sample name and condition based on parsed sample name. If sample name not found in sample map, throws error.
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
def prep_h5_conversion(h5_meta){
    // parse metadata for h5 files with matching sample name, an returns tuple of meta, h5_raw path, h5_filtered path by sample id
    def h5_data = [:]
    h5_meta.each { meta -> 
        h5_data[meta.sample][meta.input_type] =  meta.name
    }
    def h5_to_convert = []
    h5_data.each { sample, input_type ->
        h5_to_convert << [h5_meta, sample, input_type.h5_raw, input_type.h5_filtered]
    }
    return h5_to_convert
}

def parse_input_dir_src(dir_channel, src_channel, sample_map, pattern_map){
    // returns src, sample, condition, dir
    // takes list of dir paths list of dir generations sources (like cellranger, matrix (soupX), doubletFinder), sample-condition map, and parses them to create a
    // desired tuple. pattern map used to parse sample name from dir name if coming from doubletFinder or soupX, otherwise assumed to be cell ranger output and parsed as-is. sample map used to assign desired sample name and condition based on parsed sample name. If sample name not found in sample map, throws error.
    // previous runs of doubletFinder and soupX have <sample_name>_<doubletfinder/soupx> in the dir name, pattern map used for these cases, then attempt to parse sample name
    // if unable to, throw error
    // If pattern_map not matched at all, assumed cell ranger output, which stays as-is
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
        input_sample_sheet
    main:
    // group by sample, type, also strip keys with empty data
    input_by_type = input_sample_sheet.splitCsv(header: true, sep: "\t").map { metadata -> tuple (metadata.input_type, metadata.findAll{keys -> keys.value})}.view()
    input_meta_tar = input_by_type.filter { input_type, _meta -> input_type.startsWith("tar") }.map {_input_type, meta -> [meta, meta.name] }
    UNTAR_CR(
        input_meta_tar
    )
    // If there are h5 inputs, convert to CR style matrix dirs, then add to cellranger dirs
    h5_to_convert = prep_h5_conversion(input_by_type.filter { input_type, _meta -> input_type.startsWith("h5") }.map {_input_type, meta -> meta }).view()
    // h5_as_cr_dir = CONVERT_H5(h5_to_convert)
    // // add converted h5 files to input_dir_list, and set corresponding src as cellranger for parsing purposes downstream since they've converted to cellranger style output
    // input_dir_list = input_dir_list.concat(h5_as_cr_dir).collect()
    // input_dir_src_list = input_dir_src_list.concat(h5_as_cr_dir.map{ _dir -> "cellranger" }).collect()
    // parse_input_dir_src(input_dir_list, input_dir_src_list, sample_condition_map, dirname_pattern)
    // .concat(process_untar_outputs(UNTAR_CR.out, sample_condition_map, dirname_pattern))
    // .branch{ parsed_input -> 
    //         doubletfinder: parsed_input[0].toLowerCase() == "doubletfinder"
    //         matrix: parsed_input[0].toLowerCase() == "matrix"
    //         cellranger: parsed_input[0].toLowerCase() == "cellranger"
    //     }.set{src_sample_dir}
    // src_sample_dir.cellranger.view()
    // emit:
    //     doubletfinder = src_sample_dir.doubletfinder
    //     matrix = src_sample_dir.matrix 
    //     cellranger = src_sample_dir.cellranger
}