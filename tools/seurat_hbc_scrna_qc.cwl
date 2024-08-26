cwlVersion: v1.2
class: CommandLineTool
id: seurat-hbc-scrna-qc
doc: "Run custom QC on 10X output"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "pgc-images.sbgenomics.com/d3b-bixu/hbc_scrna_qc:scrnaseq"
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: seurat_hbc_scrna_qc.R
        entry:
          $include: ../scripts/seurat_hbc_scrna_qc.R
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 16000

baseCommand: [Rscript, seurat_hbc_scrna_qc.R]

inputs:
  bc_matrix_dir: { type: 'Directory?', loadListing: deep_listing, doc: "Path to 10X counts, like 'outs/raw_feature_bc_matrix'. Leave blank if using H5 inputs",
    inputBinding: { position: 1, prefix: "--data_dir"} }
  h5_matrix_inputs: { type: 'File?', doc: "Path to h5 file with 10X count results. Leave blank if using counts directory",
    inputBinding: { position: 1, prefix: "--h5_counts" } }
  sample_id: { type: string, inputBinding: { position: 1, prefix: "--sample_id"} }
  output_basename: { type: string, inputBinding: { position: 1, prefix: "--output_basename"} }
  min_umi: { type: 'int?', doc: "minimum number of umi for cell-level filtering", default: 500,
    inputBinding: { position: 1, prefix: "--min_umi"}}
  min_genes: { type: 'int?', doc: "minimum number of genes for cell-level filtering", default: 250,
    inputBinding: { position: 1, prefix: "--min_genes"} }
  min_complexity: { type: 'float?', doc: "minimum novelty score (log10GenesPerUMI)", default: 0.8,
    inputBinding: { position: 1, prefix: "--min_complexity" } }
  max_mito_ratio: { type: 'float?', doc: "maximum ratio mitochondrial reads per cell", default: 0.2,
    inputBinding: { position: 1, prefix: "--max_mito_ratio"} }
  min_gene_prevalence: { type: 'int?', doc: "Minimum number of cells a gene must be expressed in to keep after filtering", default: 10,
    inputBinding: { position: 1, prefix: "--min_gene_prevalence" } }

outputs:
  qc_barcode_metrics: { type: File, outputBinding: { glob: "*.barcode_qc.metrics.tsv"}, doc: "Table with barcodes and calculated QC metrics" }
  qc_plots: { type: File, outputBinding: { glob: "*.QC_plots.pdf"}, doc: "Pre and post filtering metrics pdf plots"}
  qc_boxplot_stats: { type: File, outputBinding: { glob: "*_summary_metrics.tsv"}, doc: "Pre and post filtering boxplot stats tsv" }
  qc_cell_counts: { type: File, outputBinding: { glob: "*.cell_counts.tsv"}, doc: "Pre and post filtering cell counts tsv" }
  qc_filtered_ct_matrix: { type: File, outputBinding: { glob: "*qc_filtered.counts_matrix.h5"}, doc: "h5 formatted counts matrix after applying minimum QC filtering" }
  qc_variable_features_plot: { type: File, outputBinding: { glob: "*.variable_features.pdf"} , doc: "PDF with a dot plot of variable genes with top 20 labeled" }
