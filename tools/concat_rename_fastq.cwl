cwlVersion: v1.2
class: CommandLineTool
id: concat-rename-fastq
doc: "Concatenate demultiplexed fastqs and rename for cell ranger input. Core count set higher than needed to speed up file download"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'ubuntu:22.04'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 16
    ramMin: 16000


baseCommand: [mkdir, -p]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.sample_name) &&
      cat
  - position: 2
    shellQuote: false
    valueFrom: >-
      > $(inputs.sample_name)/$(inputs.corrected_read_1_name) &&
      cat
  - position: 3
    shellQuote: false
    valueFrom: >-
      > $(inputs.sample_name)/$(inputs.corrected_read_2_name)


inputs:
  r1_fastqs: { type: 'File[]', doc: "Read 1 fastqs to concat",
    inputBinding: { position: 1 } }
  r2_fastqs: {type: 'File[]', doc: "Read 2 fastqs to concat",
    inputBinding: { position: 2 } }
  sample_name: {type: string, doc: "sample name, used as prefix for finding fastqs to analyze"}
  corrected_read_1_name: {type: string, doc: "Corrected read 1 name. Should be in format {Sample_name}_S1_L{3 digit lane}_R1_001"}
  corrected_read_2_name: {type: string, doc: "Corrected read 2 name. Should be in format {Sample_name}_S1_L{3 digit lane}_R2_001"}

outputs:
  renamed_dir:
    type: Directory
    outputBinding:
      glob: $(inputs.sample_name)
    doc: "Directory containing renamed fastq file"