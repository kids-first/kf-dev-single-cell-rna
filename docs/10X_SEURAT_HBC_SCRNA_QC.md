# Seurat HBC scRNA QC Outputs
Based on [this tutorial](https://github.com/hbctraining/scRNA-seq_online/blob/scRNAseq/lessons/04_SC_quality_control.md) and some Kids First-specific tweaks, the [QC step/script](../scripts/seurat_hbc_scrna_qc.R) does the following:
1. Calculate novelty score: log<sub>10</sub>(nFeaureRNA)/log<sub>10</sub>(nCount_RNA). This is described well [here](https://github.com/hbctraining/scRNA-seq_online/blob/scRNAseq/lessons/04_SC_quality_control.md#complexity)
1. Calculate mitochondrial ratio. This is done by searching gene symbols starting with `MT-` 
1. Output tsv file with QC metrics pre-filtering. Example excerpt:
   ```tsv
   cell_ids	orig.ident	nCount_RNA	nFeature_RNA	log10GenesPerUMI	mitoRatio	sample
   AAACCCAAGCAGATAT	7316-2146_SOLO_NEW_QC_OUT.EM.raw.h5	8952.9999803	4235	0.917733361994601	0.000223388808712248	7316-2146
   AAACCCAAGGATCACG	7316-2146_SOLO_NEW_QC_OUT.EM.raw.h5	1633.9999964	1122	0.949191911584887	0.000611995105387504	7316-2146
   AAACCCAAGGGTTTCT	7316-2146_SOLO_NEW_QC_OUT.EM.raw.h5	225.9999976	236	1.00798755408155	0	7316-2146
   AAACCCACAATCACGT	7316-2146_SOLO_NEW_QC_OUT.EM.raw.h5	1289.0000016	976	0.961159718063171	0	7316-2146
   AAACCCACACTACCGG	7316-2146_SOLO_NEW_QC_OUT.EM.raw.h5	1994.000028	1349	0.948567470874349	0	7316-2146
   ```

1. Apply cell-level filters. These are adjustable, and are considered a reasonable start without prior experimental knowledge. Defaults are:
   - `min_umi`: 500
   - `min_genes`: 250
   - `min_complexity`: 0.8
   - `max_mito_ratio`: 0.2
1. Drop genes with `0` counts and those present in less than `min_gene_prevalence` (default 10)
1. _Save QC-filtered counts matrix in h5 format_
1. Repeat novelty score and mito ratio calculations as well as generate post-filter QC plots.
1. Normalize read counts using `NormalizeData` and `ScaleCounts`
1. Find top 2000 variable features (transcripts) and create a dot plot labeling the top 15
1. Create pre-filter and post-filter QC plots
   - Kernel density plot of Number of transcripts (`nCount_RNA`) per cell
   - Kernel density plot of Number of genes (`nFeaureRNA`) per cell
   - Kernel density plot of  cell complexity using novelty score
   - Kernel density plot of mitochondrial counts ratio
   - Dot plot filtering effects. The intercept lines represents cell data points that will be cut out for `nUMIs` (set by `min_umi` param) and `nGenes` (set by `min_genes`) param
1. Generate pre and post-filter box plot stats and cell count stats