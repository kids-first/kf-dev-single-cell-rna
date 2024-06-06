# Seurat HBC scRNA QC Outputs
Based on [this tutorial](https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/04_SC_quality_control.md) and some Kids First-specific tweaks, the [QC step/script](../scripts/seurat_hbc_scrna_qc.R) does the following:
1. Calculate novelty score: log<sub>10</sub>(nFeaureRNA)/log<sub>10</sub>(nCount_RNA). This is described well [here](https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/04_SC_quality_control.md#complexity)
1. Calculate mitochondrial ratio. This is done by searching gene symbols starting with `MT-` 
1. _Save Seurat object with counts and metrics_
1. Create pre-filter QC plots
   - Kernel density plot of Number of transcripts (`nCount_RNA`) per cell
   - Kernel density plot of Number of genes (`nFeaureRNA`) per cell
   - Kernel density plot of  cell complexity using novelty score
   - Kernel density plot of mitochondrial counts ratio
   - Dot plot filtering effects. The intercept lines represents cell data points that will be cut out for `nUMIs` (set by `min_umi` param) and `nGenes` (set by `min_genes`) param
1. Generate prefilter box plot stats
1. Apply cell-level filters. These are adjustable, and are considered a reasonable start without prior experimental knowledge. Defaults are:
   - `min_umi`: 500
   - `min_genes`: 250
   - `min_complexity`: 0.8
   - `max_mito_ratio`: 0.2
1. Drop genes with `0` counts and those present in less than `min_gene_prevalence` (default 10)
1. Repeat novelty score and mito ratio calculations as well as generate post-filter QC plots.
1. Normalize read counts using `NormalizeData` and `ScaleCounts`
1. Find top 2000 variable features (transcripts) and create a dot plot labeling the top 15
1. _Save post-filter Seurat object with counts and metrics_
1. Generate post-filter box plot stats and cell count stats