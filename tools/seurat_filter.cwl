cwlVersion: v1.2
class: CommandLineTool
id: seurat_merge
doc: >-
  Filter Seurat QC RDS file using a SoupX RDS and Scrublet CSV

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'satijalab/seurat:4.1.0'
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: seurat_filter.R
        entry:
          $include: ../scripts/seurat_filter.R

baseCommand: [Rscript, seurat_filter.R]

inputs:
  seurat_qc_rds: { type: 'File', inputBinding: { position: 2, prefix: "--seurat_qc_rds"}, doc: "RDS object produced by Seurat QC" }
  soupx_rds: { type: 'File', inputBinding: { position: 2, prefix: "--soupx_rds"}, doc: "RDS object produced by SoupX" }
  scdblfinder_tsv: { type: 'File', inputBinding: { position: 2, prefix: "--scdblfinder_tsv"}, doc: "TSV file with doublet information produced by scDblFinder" }
  output_filename: { type: 'string?', inputBinding: { position: 2, prefix: "--output_name"}, doc: "Name for filtered Seurat RDS output" }

outputs:
  filtered_rds:
    type: File
    outputBinding:
      glob: |
        $(inputs.output_name ? inputs.output_filename : "*.rds")
    doc: "Seurat RDS filtered with mRNA contamination and doublets removed"
