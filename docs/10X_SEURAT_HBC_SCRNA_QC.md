# Seurat HBC scRNA QC Outputs
Based on [this tutorial](https://github.com/hbctraining/scRNA-seq_online/blob/scRNAseq/lessons/04_SC_quality_control.md) and some Kids First-specific tweaks, the [QC step/script](../scripts/seurat_hbc_scrna_qc.R) does the following, with output files **_emphasized_**:
1. Calculate novelty score: log<sub>10</sub>(nFeaureRNA)/log<sub>10</sub>(nCount_RNA). This is described well [here](https://github.com/hbctraining/scRNA-seq_online/blob/scRNAseq/lessons/04_SC_quality_control.md#complexity)
1. Calculate mitochondrial ratio. This is done by searching gene symbols starting with `MT-` 
1. **_qc_barcode_metrics_**:  output_basename + ".barcode_qc.metrics.tsv". Output tsv file with QC metrics pre-filtering. Example excerpt:
   ```tsv
   cell_ids	seq_source	nUMI	nGene	log10GenesPerUMI	mitoRatio	sample
   AAACCCAAGGAGAGTA	delem_test_refactor.EM.raw.h5	9118.0000543	2956	0.876462736098413	0.116979632994956	pbmc
   AAACGCTTCAGCCCAG	delem_test_refactor.EM.raw.h5	6014.999944	2103	0.879235803237773	0.0904405660955398	pbmc
   AAAGAACAGACGACTG	delem_test_refactor.EM.raw.h5	4586.999945	1778	0.887588761938559	0.0671312412671067	pbmc
   AAAGAACCAATGGCAG	delem_test_refactor.EM.raw.h5	3013.0000477	1404	0.90467602389422	0.0713574499157815	pbmc
   AAAGAACGTCTGCAAT	delem_test_refactor.EM.raw.h5	7107.0004425	2119	0.863551873376804	0.0697698282154005	pbmc
   AAAGGATAGTAGACAT	delem_test_refactor.EM.raw.h5	9537.9999555	2374	0.848226390313623	0.0849234644347971	pbmc
   ```
   Note, default Suerat column names were renamed for clarity: `orig.ident -> seq_source, nCount_RNA -> nUMI, nFeature_RNA -> nGene`

1. Apply cell-level filters. These are adjustable, and are considered a reasonable start without prior experimental knowledge. Defaults are:
   - `min_umi`: 500
   - `min_genes`: 250
   - `min_complexity`: 0.8
   - `max_mito_ratio`: 0.2
1. Drop genes with `0` counts and those present in less than `min_gene_prevalence` (default 10)
1. **_qc_filtered_ct_matrix_**: output_basename + ".qc_filtered.counts_matrix.h5". Save QC-filtered counts matrix in h5 format
1. Normalize read counts using `NormalizeData` and `ScaleCounts`
1. Find top 2000 variable features (transcripts) and create a dot plot labeling the top 15
1. **_qc_variable_features_plot_**: output_basename + ".variable_features.pdf". PDF with a dot plot of variable genes with top 20 labeled
1. **_qc_plots_**: output_basename + ".QC_plots.pdf. PDF with pre-filter and post-filter QC plots
   - Kernel density plot of Number of transcripts (`nCount_RNA`) per cell
   - Kernel density plot of Number of genes (`nFeaureRNA`) per cell
   - Kernel density plot of  cell complexity using novelty score
   - Kernel density plot of mitochondrial counts ratio
   - Dot plot filtering effects. The intercept lines represents cell data points that will be cut out for `nUMIs` (set by `min_umi` param) and `nGenes` (set by `min_genes`) param
1. **_qc_boxplot_stats_**: output_basename + ".boxplot_summary_metrics.tsv". TSV file with pre and post-filter box plot stats of the following QC metrics:
   - Number of UMIs
   - Number of Genes
   - Novelty Score
   - Mitochondrial Ratio
1. **_qc_cell_counts_**: output_basename + ".cell_counts.tsv". TSV file with number of starting cells, number of cells affected by each of the QC cutoffs, and final cell count