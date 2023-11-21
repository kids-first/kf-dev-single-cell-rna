cwlVersion: v1.2
class: CommandLineTool
id: convert-to-h5
doc: "Convert counts dir from Cell Ranger or STAR Solo to h5 file"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "pgc-images.sbgenomics.com/brownm28/soupx:1.6.2"
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: convert_to_h5.R
        entry:
          $include: ../scripts/convert_to_h5.R
  - class: ResourceRequirement
    coresMin: 2
    ramMin: 4000

baseCommand: [Rscript]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      convert_to_h5.R
inputs:
  counts_matrix_dir: { type: Directory, loadListing: deep_listing, doc: "Cell Ranger-like counts dir with matrix, features, barcodes",
    inputBinding: { position: 1, prefix: "--counts_dir" } }
  sample_name: { type: string, doc: "Sample name to put in output",
    inputBinding: { position: 1, prefix: "--sample_name" } }
  output_basename: { type: string, doc: "File name prefix name for output",
    inputBinding: { position: 1, prefix: "--output_basename" } }

outputs:
  converted_h5:
    type: File
    outputBinding:
      glob: "*.h5"
    doc: "Converted matrix dir to h5 format"

