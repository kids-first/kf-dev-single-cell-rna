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
  - position: 0
    shellQuote: false
    valueFrom: >-
      -e $(inputs.sample_name)
      -b $(inputs.barcodes.path)
      -m $(inputs.repeats.path)
      -o $(inputs.output_folder)
      $(inputs.bam.path)
      $(inputs.genes.path)

inputs:
  sample_name: {type: string, doc: "sample name, used as basename for output"}
  barcodes: {type: File[]?, doc: "list of barcodes used to filter the bam file"}
  repeats: {type: File[]?, doc: ".gtf file containing intervals to mask"}
  output_folder: {type: string, doc: "output folder"}
  bam: {type: File[]?, doc: "bam file to process"}
  genes: {type: File[]?, doc: ".gtf file with genes to analyze"}

outputs:
  velocyto_out:
    type: File
    outputBinding:
      glob: $(inputs.output_folder)/$(inputs.sample_name).loom
    doc: "Loom file, contains layered data describing calculated velocity."
