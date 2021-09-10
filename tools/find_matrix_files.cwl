cwlVersion: v1.2
class: CommandLineTool
id: find_matrix_files
doc: "From the whole cellranger output directory, return the 2 matrix files."

requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/ubuntu:18.04'

baseCommand: [cp]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -r $(inputs.count_dir.path)/filtered_gene_bc_matrices . &&
      cp -r $(inputs.count_dir.path)/raw_gene_bc_matrices .

inputs:
  count_dir: {type: Directory, doc: "Directory containing raw and filtered mtx files produced by Cell Ranger count"}

outputs:
  filtered_matrix:
    type: File
    outputBinding:
      glob: filtered_gene_bc_matrices/GRCh38/matrix.mtx.gz
  raw_matrix:
    type: File
    outputBinding:
      glob: raw_gene_bc_matrices/GRCh38/matrix.mtx.gz
