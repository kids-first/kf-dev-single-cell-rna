cwlVersion: v1.2
class: CommandLineTool
id: soupx
doc: >-
  Run SoupX on 10x output to remove mRNA contamination
  https://github.com/constantAmateur/SoupX

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "pgc-images.sbgenomics.com/d3b-bixu/soupx_r:4.1.0_SoupX"
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: run_soupX.R
        entry:
          $include: ../scripts/run_soupX.R

baseCommand: [Rscript]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      run_soupX.R --raw $(inputs.raw_matrix.path) --fil $(inputs.filtered_matrix.path) --sample_name $(inputs.sample_name) --cluster $(inputs.cluster_file.path)

inputs:
  raw_matrix:
    type: File
    doc: "H5 file raw matrix produced by Cell Ranger count"
  filtered_matrix:
    type: File
    doc: "H5 file filtered matrix produced by Cell Ranger count"
  cluster_file:
    type: File
    doc: "CSV file clusters produced by Cell Ranger count"
  sample_name:
    type: string
    doc: "Name of this sample"

outputs:
  decontaminated_matrix_dir:
    type: Directory
    outputBinding:
      loadListing: deep_listing
      glob: $(inputs.sample_name)
    doc: "Directory containing decontaminated matrix"
  decontaminated_matrix_rds:
    type: File
    outputBinding:
      loadListing: deep_listing
      glob: $(inputs.sample_name).rds
    doc: "rds containing decontaminated matrix"
  rplots:
    type: File
    outputBinding:
      glob: "Rplots.pdf"
