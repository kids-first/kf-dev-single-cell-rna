# Kids First Single-entity Workflow ANalysiS Pipeline (SWANS) Data Cleanup and QC 2.1.0:
Citation: https://doi.org/10.1101/2025.05.14.654073.
This is a NextFlow implementation of the [SWANS](https://github.com/FrancoResearchLab/SWANS) for Single Cell RNAseq Analysis.
When reviewing the original [flow diagram] (https://github.com/FrancoResearchLab/SWANS/blob/v2.1.0/markdown_images/swans.drawio.png) from the authors, this covers steps A-F.
The intent is fo the researcher to choose a set of scRNA data to set initial cutoffs, run the workflow, and review QC for any needed adjustments or simply bad data.
It can handle a mix of doubleFinder, SoupX, and CellRanger inputs, passing the inputs to the appropriate steps, and merging before initial seurat object step.
Proper FLEX input support is still forthcoming. 

## INPUTS
### Required:
Must have one of, or both a `directories` or `tar balls/h5` input
- `input_sample_sheet`: Input TSV manifest one or more of doubletFinder, soupX, cell ranger in directories, tar balls, or cellranger h5 files. Required:
  - `sample_id`: Free text string for sample name of assocaited file
  - `name`: Path/name of input file to process
  - `condition`: short (one word?) descriptor of the role of this sample in the experiment
  - `input_type`: enum describing what kind of input, must be one of:
    - `h5_raw`
    - `h5_filtered`
    - `dir_cellranger`
    - `tar_cellranger_count`
    - `tar_cellranger_multi`
    - `dir_doubletFinder`
    - `tar_doubletFinder`
    - `dir_soupX`
    - `tar_soupX`

- `project`: string for project name
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
- `soupx_start`: `no_clusters`
  - `outs`: Cell ranger outputs have clustering information
  - `no_clusters`: Have clusters calculated
  - `h5`: Input matrix is in h5 format


Refer to the [Nextflow documentation](https://www.nextflow.io/docs/latest/) for advanced usage and profiles.

## OUTPUTS
### Required
- `seurat_qc_report`: HTML QC Report. Check [here](https://github.com/FrancoResearchLab/SWANS/blob/v2.1.0/README.md#qc_report) for report guide
- `seurat_qs`: Initial seurat filtered obect required for downstream anaylsis
### If starting from scratch and flags enabled
- `doubletfinder`: DoubletFinder tar ball results. Can be reused downstream if needed
- `soupX`: soupX tar ball results. Can be used downstream and if this workflow is repeated with new params.
