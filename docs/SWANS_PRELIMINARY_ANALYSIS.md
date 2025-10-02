# SWANS 2.0: Single-entity Workflow ANalysiS Pipeline
Citation: https://doi.org/10.1101/2025.05.14.654073.
This is a NextFlow implementation of the [SWANS (Single-entity Workflow ANalysiS Pipeline)](https://github.com/FrancoResearchLab/SWANS) for Single Cell RNAseq Analysis.
It is in alpha phase, but testing on two samples using a 8 core, 64GB EC2 was successful, as well as on CAVATICA. The [test json](nf/test_inputs/test.json) has a bare bones example starting with cell ranger matrices using all defaults in the nextflow.config.
It can handle a mix of doubleFinder, SoupX, and CellRanger inputs, passing the inputs to the appropriate steps, and merging before initial seurat object step.
However, mixed mode has not been tested yet, and FLEX input support is still forthcoming. 
After running the workflow, refer to https://github.com/FrancoResearchLab/SWANS/tree/main?tab=readme-ov-file#html--shiny-report on how to use the interactive report.

## INPUTS
### Required:
- `sample_condition_map_file`: tsv with input sample, condition, remap. Last column value can be blank if original sample name to be used. Original sample name should match what is part of the directory path name
- `project`: string for project name
### Required, mutually exclusive
Parameters in which one or both must be provided, but both cannot be blank
#### Input data from _directories_, if applicable
- `input_dir_list`: Input data from one or more of doubletFinder, soupX, cell ranger in directories
- `input_dir_src_list`: Source of input dirs, options: doubletFinder, soupX, cellranger
#### Input data from _tar balls_, if applicable
- `input_dir_list`: Input data from one or more of doubletFinder, soupX, cell ranger in tar balls
- `input_dir_src_list`: Source of input tar balls, options: doubletFinder, soupX, cellranger
Must have one of, or both a `directories` or `tar balls` input
### Required with defaults:
- `organism`: `human`. Other optinon is `mouse`
- `disable_doubletfinder`: `false`. Flag to skip running doubletFinder
- `mito_cutoff`:  `15`. mitchondrial cutoff threshold, percent as int
- `ribo_cutoff`: `100`. ribosomal cutoff threshold, percent as int
- `min_feature_threshold`: `200`
- `max_feature_threshold`: `3000`
- `int_components`: `30`. number of Principal Components to use in Seurat
- `dbl_mem`: `30`. default ram in GB for doubletFinder jobs
- `dbl_threads`: `2`. default threads for doubletFinder jobs
- `disable_soupx`: `false`. Set to `true` to skip SoupX.
- `soupx_start`: `outs`
  - `outs`: Cell ranger outputs thave have clustering information
  - `no_clusters`: Have clusters calculated
  - `h5`: Input matrix is in h5 format
- `aso_processes`: `4`. Number of processes for Seurat analysis.
- `aso_memory`: `30000000000` (30 GB). Memory allocation for Seurat analysis in bytes.
- `run_azimuth`: `"n"`. Set to `"y"` to enable Azimuth reference mapping.
- `run_transferdata`: `"n"`. Set to `"y"` to enable TransferData.
- `mito_regression`: `"n"`. Set to `"y"` to enable mitochondrial regression.
- `ribo_regression`: `"n"`. Set to `"y"` to enable ribosomal regression.
- `cc_regression`: `"n"`. Set to `"y"` to enable cell cycle regression.
- `num_var_features`: `3000`. Number of variable features.
- `resolution_config`: `"0.1,0.2,0.3"`. Comma-separated clustering resolutions.
- `include_tsne`: `"n"`. Set to `"y"` to include t-SNE.
- `scale_data_features`: `"variable"`. Options:
  - `all`
  - `variable`
- `split_layers_by`: `"Sample"`. Metadata column to split layers. Options:
  - `Sample`
  - `Experiment`
- `normalization_config`: `"standard,sct"`. Comma-separated normalization methods. One or more of:
  - `sct`
  - `standard`
- `integration_config`: `"cca,rpca,harmony"`. Comma-separated integration methods. On or more of:
  - `cca`
  - `rpca`
  - `harmony`
- `ref_based_integration`: `"n"`. Set to `"y"` for reference-based integration.
- `storage`: `rds`. Output format.
- `conserved_genes`: `"n"`. Set to `"y"` to include conserved genes.
- `visualization`: `dot`. Visualization type. One or more of:
  - `feature`
  - `violin`
  - `ridge`
  - `dot`
### Optional when "parent" params is `y`
- `run_transferdata`
  - `transferdata_ref_file`: path to reference Seurat object for TransferData
  - `transferdata_reduction`: which dimensional reduction in the TRANSFERDATA REF FILE to use for TransferData
  - `transferdata_annocol`: which Seurat meta.data column to use for TransferData annotation
- `run_azimuth`
  - `azimuth_ref`: Reference preset for Azimuth Choices: 
    - `adiposeref`
    - `bonemarrowref`
    - `fetusref`
    - `heartref`
    - `humancortexref`
    - `kidneyref`
    - `lungref`
    - `pancreasref`
    - `pbmcref`
    - `tonsilref`
- `cc_regression`
  - `cc_method`: `standard`, `alternative`
- `ref_based_integration`
  - `ref_samples`: Reference samples for integration.
  - `regression_file`: Path to regression file.

Refer to the [Nextflow documentation](https://www.nextflow.io/docs/latest/) for advanced usage and profiles.

## OUTPUTS
- Analysis results tar ball. See https://github.com/FrancoResearchLab/SWANS?tab=readme-ov-file#output-structure
- Tar ball outputs of doubleFinder and SoupX if flag enabled
