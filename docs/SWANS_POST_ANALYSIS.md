# Kids First Single-entity Workflow ANalysiS Pipeline (SWANS) Final Analysis 2.1.0:
Citation: https://doi.org/10.1101/2025.05.14.654073.
This is a NextFlow implementation of the [SWANS](https://github.com/FrancoResearchLab/SWANS) for Single Cell RNAseq Analysis.
When reviewing the original [flow diagram] (https://github.com/FrancoResearchLab/SWANS/blob/v2.1.0/markdown_images/swans.drawio.png) from the authors, this covers steps I-K.
Steps L and M are skipped since nextflow provides this functionality.
See the original software description starting with ["Final (Post-Annotation) Analysis"](https://github.com/FrancoResearchLab/SWANS/blob/v2.1.0/README.md#final-post-annotation-analysis), [configs/post_annotation_configs.yaml](https://github.com/FrancoResearchLab/SWANS/blob/v2.1.0/README.md#configspost_annotation_configsyaml) and th rest of the section, and [Final Report](https://github.com/FrancoResearchLab/SWANS/blob/v2.1.0/README.md#final-report) for a sense of what this workflow does.

## INPUTS
This workflow should only be attempted only after careful review out outputs from Interactive Report.
It requires detailed knowledge of the experiment and the results
### General parameters

- `r_lib_path`: "/usr/local/lib/R/site-library/"
- `project`: Should be the same as used in Data Clean Up and QC, and Interactive Report
- `augment_prelim_config`: Provide this file as prelim_config in tar ball may be missing some entries. i.e.:
  ```
  RUN_SOUPX: true
  RUN_DOUBLETFINDER: true
  MIN_FEATURE_THRESHOLD: 200
  MAX_FEATURE_THRESHOLD: 3000
  ```

### Info about input Seurat object parameters

- `existing_analysis_tar`: Tar ball of analysis from completed interactive report run
- `user_analyzed_seurat_object_meta_sample`: 
  - default: `Sample`
  - Meta data character string to access 'Sample' in Seurat object (e.g., Samples, sample_name)
- `user_analyzed_seurat_object_meta_experiment`:
  - default: `Experiment`
  - Meta data character string to access 'Experiment' in Seurat object (e.g., Experiment, Conditions)
- `user_analyzed_seurat_object_meta_annotation`: Meta data character string holding annotation information in Seurat object (e.g., celltypes, annotation_layer)
- `user_umap_reduction`: UMAP reduction in supplied Seurat object (e.g., standard.cca.umap)
- `user_tnse_reduction`: t-SNE reduction in supplied Seurat object (e.g., standard.rpca.tsne)
- `organism`:
  default: `human`
- `annotate_provided_final_seurat_object`:
  - default: `n`
  - Annotate Seurat object? (y/n)
- `sample_list`: Sample list file from data cleanup and QC. Can be recreated following [these guidelines](https://github.com/FrancoResearchLab/SWANS/blob/v2.1.0/README.md#samplessample_list) if lost

### Final analysis parameters to apply

- `min_pct`:
  - default: `0.1`
  - Minimum percentage of cells in one of compared groups for differential gene expression
- `avg_log2fc_threshold`:
  - default: `1.5`
  - Filtering DEGs (avg_log2FC) threshold
- `final_filtering_threshold`:
  - default: `0.10`
  - Filtering threshold for adjusted p-value in pathway/GSEA results
- `final_seurat_normalization_method`:
  - default: `standard`
  - Final normalization method for Seurat (standard, sct)
- `final_seurat_integration_method`:
    - default: `harmony`
    - Final integration method for Seurat (cca, rpca, harmony)
- `final_resolution`:
    - default: `0.1`
    - Final resolution value for clustering
- `cluster_annotation_file`: TSV file with cluster names and cell type assignments
- `final_storage`:
    - default: `rds`
    - Final storage format (rds, sceo, cloupe, cellchat, cellphone)
- `final_user_gene_file`: Relative path with filename to file containing genes of interest
- `final_visualization`:
    - default: `dot`
    - Visualization type for final_user_gene_file (feature, violin, ridge, dot)
- `final_conserved_genes`:
    - default: `n`
    - Run conserved genes analysis? (y/n) - Only run with multiple experiments
- `final_threads`:
    - default: `30`
    - Number of threads for computation
- `memory_mb`:
    - default: `30000`
    - Memory in MB for computation

### Run trajectory analysis
- `run_trajectory_analysis`:
  - default: `false`
  - Run Monocle3 trajectory analysis? (boolean)
- `partition_trajectory`:
  - default: `y`
  - If `y`, CLUSTER_ANNOTATION_FILE needs a 'partition' column with a value for each cluster

