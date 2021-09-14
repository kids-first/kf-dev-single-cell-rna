cwlVersion: v1.2
class: CommandLineTool
id: add_filtered_matrices
doc: "Add the doublet filtered matrices to the cellranger output directories."

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
      -r $(inputs.count_dir.path) . &&
      cp $(inputs.raw_matrix.path) $(inputs.count_dir.basename)/raw_gene_bc_matrices/GRCh38/matrix.mtx &&
      cp $(inputs.filtered_matrix.path) $(inputs.count_dir.basename)/filtered_gene_bc_matrices/GRCh38/matrix.mtx &&
      gzip --force $(inputs.count_dir.basename)/raw_gene_bc_matrices/GRCh38/matrix.mtx &&
      gzip --force $(inputs.count_dir.basename)/filtered_gene_bc_matrices/GRCh38/matrix.mtx

inputs:
  count_dir: {type: Directory, doc: "Directory containing raw and filtered mtx files produced by Cell Ranger count"}
  filtered_matrix: {type: File, doc: "doublet removed filtered matrix"}
  raw_matrix: {type: File, doc: "doublet removed raw matrix"}

outputs:
  out_dir:
    type: Directory
    outputBinding:
      glob: $(inputs.count_dir.basename)
