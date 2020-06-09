cwlVersion: v1.0
class: CommandLineTool

id: velocyto_run

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/velocyto:0.17.17'
  - class: InlineJavascriptRequirement

baseCommand: [velocyto, run]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -e $(inputs.sample_name)
      -b $(inputs.barcodes)
      -m $(inputs.repeats)
      -o $(inputs.output_folder)
      $(inputs.bam)
      $(inputs.genes)

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
      glob: *.loom
