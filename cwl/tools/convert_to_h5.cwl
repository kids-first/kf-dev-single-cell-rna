cwlVersion: v1.2
class: CommandLineTool
id: convert-to-h5
doc: "Convert counts dir from STAR Solo to h5 file. If user wishes to use a mutli-map count raw file from the results, use raw_count_choice"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "pgc-images.sbgenomics.com/d3b-bixu/soupx_r:1.6.2"
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing: 
      - $(inputs.solo_counts_dir)
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
      mkdir RAW
      && cp $(inputs.solo_counts_dir.path)/$(inputs.soloFeatures)/raw/*.tsv RAW/ 
      && cp $(inputs.solo_counts_dir.path)/$(inputs.soloFeatures)/raw/$(["Unique",null].indexOf(inputs.raw_count_choice) != -1 ? "matrix" : "UniqueAndMult-" + inputs.raw_count_choice).mtx RAW/matrix.mtx &&
  - position: 1
    shellQuote: false
    valueFrom: >-
      Rscript convert_to_h5.R
      --counts_dir "RAW/"
      --sample_name $(inputs.sample_name)
      --output_basename $(inputs.output_basename).$(inputs.raw_count_choice).raw
  - position: 2
    shellQuote: false
    valueFrom: >-
      && Rscript convert_to_h5.R
      --counts_dir $(inputs.solo_counts_dir.path)/$(inputs.soloFeatures)/filtered/
      --sample_name $(inputs.sample_name)
      --output_basename $(inputs.output_basename).filtered

inputs:
  solo_counts_dir: { type: Directory, loadListing: deep_listing, doc: "Cell Ranger-like counts dir with matrix, features, barcodes" }
  sample_name: { type: string, doc: "Sample name to put in output" }
  soloFeatures: { type: [ 'null', {type: enum, name: soloFeatures, symbols: ["Gene", "SJ", "GeneFull", "GeneFull_ExonOverIntron", "GeneFull_Ex50pAS"]}], doc: "This opt used in STAR Solo will determine the folder structure",
    default: "GeneFull"}
  raw_count_choice: {type: ['null', {type: enum, name: raw_count_choice, symbols: ["no_mm", "Uniform", "PropUnique", "EM", "Rescue"]}], doc: "Based on `soloMultiMappers`, if you wish to include/handle multi-gene hits in downstream anaylsis instead of default (ignore multi-gene mappers), pick the method you want to use",
      default: "no_mm" }
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

