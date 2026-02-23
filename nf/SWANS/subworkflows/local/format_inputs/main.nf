#!/usr/bin/env nextflow

include { UNTAR_CR } from '../../../modules/local/untar/main.nf'
include { CONVERT_H5 } from '../../../modules/local/convert_h5/main.nf'


workflow format_inputs {
    take:
        input_sample_sheet
    main:
    // Convert sample_sheet to queue channel of tuples(input_type, meta)
    input_by_type = input_sample_sheet.splitCsv(header: true, sep: "\t").map { metadata -> tuple(metadata.input_type, metadata.findAll{keys -> keys.value})}
    input_meta_tar = input_by_type.filter { input_type, _meta -> input_type.startsWith("tar") }
        .map {_input_type, meta -> [meta, meta.name] }
        .unique { _meta, filename -> filename }
    UNTAR_CR(
        input_meta_tar.map{ meta, file -> [meta + ["input_type": meta.input_type.replace("tar_", "dir_")], file] }
    )
    // Correct the metadata for the sample IDs matching tar_cellranger_multi, if any
    tar_cellranger_multi_meta = input_by_type.filter { input_type, _meta -> input_type == "tar_cellranger_multi" }
        .map {_input_type, meta -> tuple(meta.subMap('name', 'sample_id'), meta) }
    // if input_type is dir_cellranger_multi, then it's a multi, so replace with meta from tar CR multi, update meta.name
    // else process as normal, but update associated path with dir path from untar
    UNTAR_CR.out.branch { meta, dir_path ->
        multi: meta.input_type.equals("dir_cellranger_multi")
            return dir_path.map{ dir -> [["name": meta.name, "sample_id": dir.baseName], dir] }
        count: meta.input_type.equals("dir_cellranger_count")
            return [meta, dir_path]
        doublet: meta.input_type.equalsIgnoreCase("dir_doubletFinder")
            return [meta, dir_path]
        soupx: meta.input_type.equalsIgnoreCase("dir_soupX")
            return [meta, dir_path]
        }.set {per_sample_untar}
    // if there were tar metadata, prep to associate untarred multi dirs with the correct sample metadata, else initialize empty channel
    multi_dir = tar_cellranger_multi_meta
        .join(per_sample_untar.multi.flatMap())
        .map{_joinkey, meta, untar_path ->
            meta.name = untar_path
            tuple(meta, untar_path)
        }

    // If there are h5 inputs, convert to CR style matrix dirs, then add to cellranger dirs
    h5_to_convert = input_by_type.filter { input_type, _meta -> input_type.startsWith("h5") }
        .map {_input_type, meta -> [meta.sample_id, meta] }
        .groupTuple()
        .map { _sample_id, meta_list -> tuple(
            meta_list[0] + ["input_type": "dir_cellranger"], meta_list.find{ sub_meta -> sub_meta.input_type == "h5_raw" }.name, meta_list.find{ sub_meta -> sub_meta.input_type == "h5_filtered" }.name
            )
        }
    h5_as_cr_dir = CONVERT_H5(h5_to_convert)
    // start collating cellranger, soupX. and doubletFinder dirs as applicable from input dirs, tars, and h5 conversions
    cellranger = input_by_type.filter { input_type, _meta -> input_type == "dir_cellranger" }
        .map {_input_type, meta -> tuple(meta, meta.name) }
        .concat(h5_as_cr_dir)
        .concat(per_sample_untar.count)
        .concat(multi_dir)
    doubletFinder = input_by_type.filter { input_type, _meta -> input_type == "dir_doubletFinder" }
        .map {_input_type, meta -> tuple(meta, meta.name) }
        .concat(per_sample_untar.doublet)
    soupX = input_by_type.filter { input_type, _meta -> input_type == "dir_soupX" }
        .map {_input_type, meta -> tuple(meta, meta.name) }


    emit:
        cellranger
        doubletFinder
        soupX
}
