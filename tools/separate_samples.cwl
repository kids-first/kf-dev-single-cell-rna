cwlVersion: v1.2
class: CommandLineTool
id: separate_samples
doc: "Take a directory with fastq files and create array of fastq1s, fastq2s, and sample names."

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/loompy:2.0.16'
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: separate_samples.py
        entry:
          $include: ../scripts/separate_samples.py

baseCommand: [python]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      separate_samples.py --dir $(inputs.input_dir.path) > samples.txt

inputs:
  input_dir: {type: 'Directory', doc: "Directory with fastqs"}
  paired: {type: "boolean?", inputBinding: {position: 1, prefix: -p}, doc: "Flag for paired samples"}

outputs:
  fastq1s:
    type: File[]
    outputBinding:
      glob: fastq1/*.fastq.gz
    doc: "Array of fastq 1s"
  fastq2s:
    type: File[]?
    outputBinding:
      glob: fastq2/*.fastq.gz
    doc: "Array of fastq 2s"
  samples:
    type: File
    outputBinding:
      glob: samples.txt
    doc: "File with sample names"
