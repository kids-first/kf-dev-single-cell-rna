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
      seurat_merge.R --matrix_files ${var name_string = ""; for (var i=0; i<inputs.matrix_rds_files.length-1; i++){name_string += inputs.matrix_rds_files[i].path+","} name_string += inputs.matrix_rds_files[inputs.matrix_rds_files.length-1].path; return name_string} --output_name $(inputs.output_name)

inputs:
  matrix_rds_files:
    type: File[]
    doc: "RDS files containing count matrices filtered by SoupX"
  output_name:
    type: string
    doc: "Name with which to tag output matrix"

outputs:
  merged_matrix:
    type: File
    outputBinding:
      glob: $(inputs.output_name).merged.RDS
    doc: "RDS file containing merged count matrix"
