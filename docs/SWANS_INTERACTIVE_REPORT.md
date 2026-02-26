# Kids First Single-entity Workflow ANalysiS Pipeline (SWANS) Interactive Report 2.1.0:
Citation: https://doi.org/10.1101/2025.05.14.654073.
This is a NextFlow implementation of the [SWANS](https://github.com/FrancoResearchLab/SWANS) for Single Cell RNAseq Analysis.
When reviewing the original [flow diagram] (https://github.com/FrancoResearchLab/SWANS/blob/v2.1.0/markdown_images/swans.drawio.png) from the authors, this covers steps G-H.
See the original software description starting with ["clustering and DGE analysis"](https://github.com/FrancoResearchLab/SWANS/blob/v2.1.0/README.md#preliminary-analysis) for a sense of what this workflow does.

## INPUTS
### Required:
- `project`: string for project name. _Must match name used in `initial_seurat_object`_
- `initial_seurat_object`: qs file output from SWANS Data Cleanup and QC
### Required with defaults
- `organism`: organism data is from. Can be human or mouse. Default is `human`
- `mito_cutoff`: mitochondrial cutoff percent as int. Set to 100 to filter nothing. Default is `15`
- `ribo_cutoff`: ribosomal cutoff percent as int. Set to 100 to filter nothing. Default is `100`
- `int_components`: number of pcs to use in seurat. Default is `30`
- `aso_processes`: number of processes to use for parallel processing in analyze_seurat_object module. Default is `4`
- `aso_memory`: memory allocation for analyze_seurat_object. Default is `30000000000`
- `run_azimuth`: whether to run Azimuth annotation. Default is `n`
- `azimuth_ref`: azimuth reference dataset to use for annotation
  - human ref options: adiposeref, bonemarrowref, fetusref, heartref, humancortexref, kidneyref, lungref, pancreasref
  - mouse ref options:  mousecortexref
- `run_transferdata`: whether to run TransferData. Default is `n`
- `mito_regression`: whether mitochondrial percent should be regressed. Default is `n`
- `ribo_regression`: whether ribosomal percent should be regressed. Default is `n`
- `num_var_features`: number of variable features (genes) to use in normalization. Default is `3000`
- `resolution_config`: resolution value(s) as csv string. Default is `0.1,0.2,0.3`
- `include_tsne`: whether to make tsne plots. Default is `n`
- `scale_data_features`: which features to use for scale data (variable or all). Default is `variable`
- `split_layers_by`: metadata by which to split the object into layers. Default is `Sample`
- `normalization_config`: normalization method for Seurat as csv string. Default is `standard,sct`
- `integration_config`: integration method for seurat as csv. Default is `cca,rpca,harmony`
- `ref_based_integration`: whether to use reference-based integration. Default is `n`
- `storage`: storage format for output (rds or sceo). Default is `rds`
- `conserved_genes`: whether to run conserved genes analysis. Default is `n`
- `visualization`: visualization type for plots. Default is `dot`
### Optional/dependent params
Params required only if corresponding flag is `y`
- `transferdata_ref_file`: path to reference Seurat object for TransferData
- `transferdata_reduction`: which dimensional reduction in the reference file to use for transferdata
- `transferdata_annocol`: which seurat meta.data column to use for transferdata annotation
- `cc_method`: method for cell cycle scoring (standard or alternative)
- `ref_samples`: which sample(s) to use as references for ref-based integration
### Optional params
- `regression_file`: file of additional genes to be regressed
- `user_gene_file`: relative path with filename to file containing genes of interest
## Outputs
 - `analysis_result`: Tarball with all images produced for shiny report. Download, untar, and see [docs](https://github.com/FrancoResearchLab/SWANS/blob/v2.1.0/README.md#html--shiny-report) for use
