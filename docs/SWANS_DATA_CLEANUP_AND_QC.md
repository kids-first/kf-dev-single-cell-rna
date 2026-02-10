# Kids First Single-entity Workflow ANalysiS Pipeline (SWANS) Data Cleanup and QC 2.1.0:
Citation: https://doi.org/10.1101/2025.05.14.654073.
This is a NextFlow implementation of the [SWANS](https://github.com/FrancoResearchLab/SWANS) for Single Cell RNAseq Analysis.
When reviewing the original [flow diagram] (https://github.com/FrancoResearchLab/SWANS/blob/v2.1.0/markdown_images/swans.drawio.png) from the authors, this covers steps A-F.
The intent is fo the researcher to choose a set of scRNA data to set initial cutoffs, run the workflow, and review QC for any needed adjustments or simply bad data.
It can handle a mix of doubleFinder, SoupX, and CellRanger inputs, passing the inputs to the appropriate steps, and merging before initial seurat object step.
Proper FLEX input support is still forthcoming. 

## INPUTS
### Required:
- `sample_condition_map_file`: tsv with input sample, condition, remap. Last column value can be blank if original sample name to be used. Original sample name should match what is part of the directory path name
- `project`: string for project name
### Required, mutually exclusive
Parameters in which one or both must be provided, but both cannot be blank.
**Note**: `matrix` typically manes `soupX`
#### Input data from _directories_, if applicable
- `input_dir_list`: Input data from one or more of doubletFinder, matrix, cell ranger in directories
- `input_dir_src_list`: Source of input dirs, options: doubletFinder, matrix, cellranger
#### Input data from _tar balls_, if applicable
- `input_tar_list`: Input data from one or more of doubletFinder, matrix, cell ranger in tar balls
- `input_tar_src_list`: Source of input tar balls, options: doubletFinder, matrix, cellranger
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


Refer to the [Nextflow documentation](https://www.nextflow.io/docs/latest/) for advanced usage and profiles.

## OUTPUTS
### Required
- `seurat_qc_report`: HTML QC Report. Check [here](https://github.com/FrancoResearchLab/SWANS/blob/v2.1.0/README.md#qc_report) for report guide
- `seurat_qs`: Initial seurat filtered obect required for downstream anaylsis
- `initial_seurat_qc`: Tar ball of outputs from step `F`.  This is required for downstream steps of SWANS. Structured as follows:
  ```
  `project`_initial_seurat_qc
  ├── RDS
  │   └── `project`_initial_seurat_object.qs
  ├── figures
  │   ├── `project`_qc_1_scatter1.png
  │   ├── `project`_qc_1_scatter2.png
  │   ├── `project`_qc_1_vln.png
  │   ├── `project`_qc_2_scatter1.png
  │   ├── `project`_qc_2_scatter2.png
  │   └── `project`_qc_2_vln.png
  ├── report
  │   ├── `project`_samples.sample_list
  │   ├── prelim_configs.yaml
  │   └── qc_report
  │       └── `project`_qc_report.html
  └── tables
  ```
### If starting from scratch and flags enabled
- `doubletfinder`: DoubletFinder tar ball results. Can be reused downstream if needed
- `soupX`: soupX tar ball results. Can be used downstream and if this workflow is repeated with new params.
