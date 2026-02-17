#!/usr/bin/env nextflow

include { UNTAR_CR } from '../../../modules/local/untar/main.nf'
include { CONVERT_H5 } from '../../../modules/local/convert_h5/main.nf'


workflow format_inputs {
    take:
        input_sample_sheet
    main:
    // group by sample, type, also strip keys with empty data, use first meta for tar flex, add meta later
    input_by_type = input_sample_sheet.splitCsv(header: true, sep: "\t").map { metadata -> tuple (metadata.input_type, metadata.findAll{keys -> keys.value})}
    input_meta_tar = input_by_type.filter { input_type, _meta -> input_type.startsWith("tar") }
        .map {_input_type, meta -> [meta, meta.name] }
        .unique { _meta, filename -> filename }
    UNTAR_CR(
        input_meta_tar
    )
    UNTAR_CR.out.view()
//     // If there are h5 inputs, convert to CR style matrix dirs, then add to cellranger dirs
//     h5_to_convert = input_by_type.filter { input_type, _meta -> input_type.startsWith("h5") }
//         .map {_input_type, meta -> [meta.sample_id, meta] }
//         .groupTuple(by: 0)
//         .map { _sample_id, meta_list -> tuple(
//             meta_list[0], meta_list.find{ sub_meta -> sub_meta.input_type == "h5_raw" }.name, meta_list.find{ sub_meta -> sub_meta.input_type == "h5_filtered" }.name
//             )
//         }
//     h5_as_cr_dir = CONVERT_H5(h5_to_convert)
//     // start collating cellranger, soupX. and doubletFinder dirs as applicable from input dirs, tars, and h5 conversions
//     cr_from_tar = UNTAR_CR.out.filter { meta, _parsed_id, _dir_path -> meta.input_type == "dir_cellranger" }
//         .map {meta, _parsed_id, dir_path -> tuple(meta, dir_path) }
//     cellranger = input_by_type.filter { input_type, _meta -> input_type == "dir_cellranger" }
//         .map {_input_type, meta -> tuple(meta, meta.name) }
//         .concat(h5_as_cr_dir)
//         .concat(cr_from_tar)


//     emit:
//         // doubletfinder = src_sample_dir.doubletfinder
//         // matrix = src_sample_dir.matrix
//         cellranger
}