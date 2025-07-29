#!/usr/bin/env nextflow

include { DOUBLETFINDER } from './modules/local/doubletFinder/main.nf'
include { SOUPX } from './modules/local/soupX/main.nf'
include { COLLATE_OUTPUTS } from './modules/local/collate_outputs/main.nf'
include { CREATE_INITIAL_SEURAT } from './modules/local/create_initial_seurat/main.nf'

workflow {
    main:
    sample_list = Channel.fromList(params.sample_list)
    condition_list = Channel.fromList(params.condition_list).collect()
    starting_data = Channel.value(params.starting_data)
    input_dir_list = Channel.fromPath(params.input_dir_list, type: 'dir')
    data_dir = "data/endpoints"

    input_list = sample_list.merge(input_dir_list).map { sample, input_dir -> [sample, input_dir]}
    if (!params.disable_doubletfinder){
        DOUBLETFINDER(
            input_list,
            starting_data,
            data_dir,
            params.mito_cutoff,
            params.min_feature_threshold,
            params.int_components,
            params.organism
        )
    }
    if (!params.disable_soupx){
        SOUPX(
            input_list,
            params.data_type,
            data_dir,
            starting_data
        )
    }
    // collate results so that next step can use them all together
    samples = SOUPX.out.map { it[0] }.collect().combine(DOUBLETFINDER.out.map { it[0] }.collect())
    input_dirs = SOUPX.out.map { it[1] }.collect().combine(DOUBLETFINDER.out.map { it[1] }.collect())
    COLLATE_OUTPUTS(
        samples,
        input_dirs,
        data_dir
    )
    COLLATE_OUTPUTS.out.view()
    seurat_filename = "data/endpoints/$params.project/analysis/RDS/${params.project}_initial_seurat_object.qs"
    CREATE_INITIAL_SEURAT(
        COLLATE_OUTPUTS.out,
        sample_list.collect(),
        condition_list,
        input_dir_list.collect(),
        params.disable_soupx ? "cellranger" : "soupX",
        'y',
        params.mito_cutoff,
        params.ribo_cutoff,
        params.min_feature_threshold,
        params.max_feature_threshold,
        seurat_filename
    )

}