cwlVersion: v1.2
class: CommandLineTool
id: star-solo-mimic-cr
doc: "Convert dir output of STAR Solo and add in clusters file to mimic cell ranger count output format"

requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 2
    ramMin: 4000

baseCommand: [mkdir, -p]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      $(inputs.output_name)/outs/raw_feature_bc_matrix
      && mkdir $(inputs.output_name)/outs/filtered_feature_bc_matrix
      && cp $(inputs.solo_counts_dir.path)/raw/* $(inputs.output_name)/outs/raw_feature_bc_matrix/
      && cp $(inputs.solo_counts_dir.path)/filtered/* $(inputs.output_name)/outs/filtered_feature_bc_matrix/
      && gzip $(inputs.solo_counts_dir.path)/raw/* $(inputs.output_name)/outs/raw_feature_bc_matrix/*
      && gzip $(inputs.solo_counts_dir.path)/filtered/* $(inputs.output_name)/outs/filtered_feature_bc_matrix/*
inputs:
  solo_counts_dir: { type: Directory, loadListing: deep_listing, doc: "Subdirectory, like 'GeneFull' from STAR Solo with raw and filtered counts matrices" }
  output_name: { type: string, doc: "Directory name for output" }

outputs:
  cr_like_counts_dir:
    type: Directory
    outputBinding:
      glob: "$(inputs.output_name)"
    doc: "Resultant CR-like directory"
