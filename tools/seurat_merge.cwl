cwlVersion: v1.2
class: CommandLineTool
id: seurat_merge
doc: >-
  Merge multiple 10X single-cell count matrices using Seurat

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/soupx_r:4.1.0_SoupX'
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: seurat_merge.R
        entry:
          $include: ../scripts/seurat_merge.R

baseCommand: [Rscript]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      seurat_merge.R --matrix_dirs $(inputs.matrix_rds_files.map(function(i){return i.path}).join()) \
      --output_name $(inputs.output_name)

inputs:
  matrix_dirs:
    type: Directory[]
    doc: "Directories containing count matrices filtered by SoupX"
  doublets_files:
    type: File[]
    doc: "Csv files with barcodes, doublet score, and predicted doublet boolean"
  output_name:
    type: string
    doc: "Name with which to tag output matrix"

outputs:
  merged_matrix:
    type: File
    outputBinding:
      glob: $(inputs.output_name).merged_matrix.RDS
    doc: "RDS file containing merged count matrix"
  merged_object:
    type: File
    outputBinding:
      glob: $(inputs.output_name).merged_object.RDS
    doc: "RDS file containing merged Seurat object"
