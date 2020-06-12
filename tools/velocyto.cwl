cwlVersion: v1.0
class: CommandLineTool
id: velocyto_run
doc: "Run velocyto to compute RNA velocity of a single cell RNA dataset."

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/velocyto:0.17.17'
  - class: ResourceRequirement
    ramMin: 40000
  - class: InlineJavascriptRequirement

baseCommand: [velocyto, run]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -e $(inputs.sample_name)
  - position: 2
    shellQuote: false
    valueFrom: >-
      -b $(inputs.barcodes.path)
  - position: 3
    shellQuote: false
    valueFrom: >-
      -m $(inputs.repeats.path)
  - position: 4
    shellQuote: false
    valueFrom: >-
      -o $(inputs.output_folder)
  - position: 5
    shellQuote: false
    valueFrom: >-
      $(inputs.bam.path)
  - position: 6
    shellQuote: false
    valueFrom: >-
      $(inputs.genes.path)

inputs:
  sample_name:
    type: string
  barcodes:
    type: File
  repeats:
    type: File
  output_folder:
    type: string
  bam:
    type: File
  genes:
    type: File

outputs:
  velocyto_out:
    type: File
    outputBinding:
      glob: $(inputs.output_folder)/$(inputs.sample_name).loom
