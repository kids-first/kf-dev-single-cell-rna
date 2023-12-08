cwlVersion: v1.2
class: CommandLineTool
id: convert-to-h5
doc: "Convert counts dir from STAR Solo to h5 file"

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

baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      Rscript convert_to_h5.R
      --counts_dir $(inputs.solo_counts_dir.path)/$(inputs.soloFeatures)/filtered/
      --sample_name $(inputs.sample_name)
      --output_basename $(inputs.output_basename).filtered
  - position: 1
    shellQuote: false
    valueFrom: >-
      && Rscript convert_to_h5.R
      --counts_dir $(inputs.solo_counts_dir.path)/$(inputs.soloFeatures)/raw/
      --sample_name $(inputs.sample_name)
      --output_basename $(inputs.output_basename).raw

inputs:
  solo_counts_dir: { type: Directory, loadListing: deep_listing, doc: "Cell Ranger-like counts dir with matrix, features, barcodes" }
  sample_name: { type: string, doc: "Sample name to put in output" }
  soloFeatures: { type: [ 'null', {type: enum, name: soloFeatures, symbols: ["Gene", "SJ", "GeneFull", "GeneFull_ExonOverIntron", "GeneFull_Ex50pAS"]}], doc: "This opt used in STAR Solo will determine the folder structure",
    default: "GeneFull"}
  output_basename: { type: string, doc: "File name prefix name for output" }

outputs:
  raw_converted_h5:
    type: File
    outputBinding:
      glob: "*raw.h5"
    doc: "Converted raw matrix dir to h5 format"
  filtered_converted_h5:
    type: File
    outputBinding:
      glob: "*filtered.h5"
    doc: "Converted filtered matrix dir to h5 format"

