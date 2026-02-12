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
def prep_h5_conversion(input_file_src_list, input_file_list, sample_map){
    // takes list of input file sources and list of input files, filters for h5 files, and creates channel of src, sample, condition, file for h5 inputs to be converted to count matrices. Sample and condition parsed from sample map based on file name. If file name not found in sample map, throws error.
    return input_file_src_list.merge(input_file_list) { src, file -> [src, file] }
    .filter { src, _file -> src.startsWith("h5") }
    .map { src, file ->
        def sname = file.name.replaceAll(/\.cellranger.*\.h5$/, "")
        if (sample_map.containsKey(sname)){
            def sample_name = sample_map[sname][1] ?: sname
            return [sample_name, src, file]
        } else {
            error("Sample name ${sname} not found in sample_condition_map_file. Please ensure all samples are mapped.")
        }
    }
}

def parse_input_dir_src(dir_channel, src_channel, sample_map, pattern_map){
    // takes list of dir paths list of dir generations sources (like cellranger, matrix (soupX), doubletFinder), sample-condition map, and parses them to create a
    // unified channel of src, sample, condition, dir. pattern map used to parse sample name from dir name if coming from doubletFinder or soupX, otherwise assumed to be cell ranger output and parsed as-is. sample map used to assign desired sample name and condition based on parsed sample name. If sample name not found in sample map, throws error.
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
        input_file_src_list
        input_file_list
        sample_condition_map_file
        input_dir_list
        input_dir_src_list
    main:
    // dir names typically drive sample names, but not always desired. use sample map to enforce desired names
    // filter out h5 files
    input_meta_tar = input_file_src_list.merge(input_file_list) { src, tar -> [src, tar] }.filter { src, _tar -> !src.startsWith("h5") }
    input_meta_tar.view()
    UNTAR_CR(
        input_meta_tar
    )
    // parse sample map file
    sample_condition_map =  parse_map_file(sample_condition_map_file.text)
    // initialize dir, tar, and src data as one channel with src, sample_name, condition, dir
    dirname_pattern = [
        doubletfinder: /([^\/]+)[_\/]doubletFinder/,
        matrix: /([^\/]+)[_\/]soupX/
    ]
    // If there are h5 inputs, convert to CR style matrix dirs, then add to cellranger dirs
    h5_to_convert = prep_h5_conversion(input_file_src_list, input_file_list, sample_condition_map).groupTuple(by: 0)
    h5_as_cr_dir = CONVERT_H5(h5_to_convert)
    // add converted h5 files to input_dir_list, and set corresponding src as cellranger for parsing purposes downstream since they've converted to cellranger style output
    input_dir_list = input_dir_list.concat(h5_as_cr_dir).collect()
    input_dir_list.view()
    input_dir_src_list = input_dir_src_list.concat(h5_as_cr_dir.map{ _dir -> "cellranger" }).collect()
    input_dir_src_list.view()
    parse_input_dir_src(input_dir_list, input_dir_src_list, sample_condition_map, dirname_pattern)
    .concat(process_untar_outputs(UNTAR_CR.out, sample_condition_map, dirname_pattern))
    .branch{ parsed_input -> 
            doubletfinder: parsed_input[0].toLowerCase() == "doubletfinder"
            matrix: parsed_input[0].toLowerCase() == "matrix"
            cellranger: parsed_input[0].toLowerCase() == "cellranger"
        }.set{src_sample_dir}
    emit:
        doubletfinder = src_sample_dir.doubletfinder
        matrix = src_sample_dir.matrix 
        cellranger = src_sample_dir.cellranger
}