# 10X Refinement Workflow

The workflow script that runs the tools is `workflows/kf_single_cell_10x_refinement.cwl`

[SoupX](https://github.com/constantAmateur/SoupX) is used for subtraction of the RNA background.
[Scrublet](https://github.com/swolock/scrublet) is used to score and predict doublets.
Decontaminated outputs are aggregated using the [Seurat](https://satijalab.org/seurat/) R package from the Satija lab at the New York Genome Center.
Original workflow design heavily contributed to by Erin Reichenbee of DBHi.

## Software

- soupX 4.1.0
- scrublet 0.2.3
- Seurat 4.0.4

## Inputs
### multi-step
 - `output_basename`: basename used to name output files
 - `sample_name`: used as prefix for finding fastqs to analyze, e.g. 1k_PBMCs_TotalSeq_B_3p_LT_antibody if the names of the underlying fastqs are of the form 1k_PBMCs_TotalSeq_B_3p_LT_antibody_S1_L001_I1_001.fastq.gz, one per input fastq in the same order
### soupX
 - `cellranger_matrix_raw`: Raw feature matrix file from Cellranger
 - `cellranger_matrix_filtered`: Filtered feature matrix file from Cellranger
 - `cellranger_cluster`: CSV containing cluster information from Cellranger
### scrublet
 - `expected_doublet_rate`: expected doublet rate, usually specific to the method; default 0.06 for 10X
 - `doublet_score_threshold`: doublet cut-off, cells with greater scores will be labelled as doublets; must be between 0 and 1
 - `count_min`: minimum expression count to retain a gene
 - `cell_min`: minimum number of cells a gene must be in to be retained
 - `min_gene_variability_pctl`: Keep the most highly variable genes (in the top min_gene_variability_pctl percentile), as measured by the v-statistic
 - `n_prin_comps`: Number of PCs to use for clustering
 - `ram`: In GB
 - `cpus`: Number of CPUs to request

### Outputs
- `soupx_rplots`: PDF R plot made by soupX 
- `scrublet_histogram`: PNG histogram made by scrublet 
- `scrublet_doublets`: CSV containing expression matrix with doublets removed by scrublet
- `decontam_matrix`: RDS file containing merged count matrix from Seurat 
- `decontam_object`: RDS file containing merged Seurat R object
