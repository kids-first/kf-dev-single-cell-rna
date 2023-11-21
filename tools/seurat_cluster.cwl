cwlVersion: v1.2
class: CommandLineTool
id: seurat-cluster
doc: "Convert counts dir from Cell Ranger or STAR Solo to h5 file"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "pgc-images.sbgenomics.com/brownm28/soupx:1.6.2"
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: seurat_cluster.R
        entry:
          $include: ../scripts/seurat_cluster.R
  - class: ResourceRequirement
    coresMin: 2
    ramMin: 4000

baseCommand: [Rscript]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      seurat_cluster.R
inputs:
  filtered_counts_matrix_dir: { type: 'Directory?', loadListing: deep_listing, doc: "Cell Ranger-like raw counts dir with matrix, features, barcodes if desired over h5",
    inputBinding: { position: 1, prefix: "--data" } }
  filtered_counts_h5: { type: 'File?', doc: "Cell Ranger-like raw counts h5, if desired over matrix dir",
    inputBinding: { position: 1, prefix: "--data" } }
  sample_name: { type: string, doc: "Sample name to put in output",
    inputBinding: { position: 1, prefix: "--name" } }
  min_features: { type: 'int?', doc: "Minimum number of genes observed in a cell to retain", default: 200,
    inputBinding: { position: 1, prefix: "--min_features"} }
  retain_features: { type: 'int?', doc: "Number of most-variable features to initially retain", default: 2000,
    inputBinding: { position: 1, prefix: "--retain_features"} }
  norm_method: { type: ['null', {type: enum, name: norm_method, symbols: ["LogNormalize", "CLR", "RC"]}],
    default: "LogNormalize", doc: "Normalization to apply to counts LogNormalize, CLR, RC" }
  num_pcs: { type: "int?", doc: "Minimum number of principal components to retain for clustering", default: 10,
    inputBinding: { position: 1, prefix: "--num_pcs" } }
  pc_cut: { type: "float?", doc: "p-value cutoff for determing the number of prinicipal components to retain for clustering", default: 0.05,
    inputBinding: { position: 1, prefix: "--pc_cut" } }
  knn_granularity: { type: "float?", doc: "KNN clustering granularity parameter", default: 0.5,
    inputBinding: { position: 1, prefix: "--knn_granularity" } }
  output_basename: { type: string, doc: "File name prefix name for output",
    inputBinding: { position: 1, prefix: "--output_basename" } }

outputs:
  clusters_csv:
    type: File
    outputBinding:
      glob: "*.csv"
    doc: "SoupX-comaptible clusters file"

