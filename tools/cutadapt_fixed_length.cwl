cwlVersion: v1.2
class: CommandLineTool
id: cutadapt-fixed-length
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/cutadapt:3.4'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 16000

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      cutadapt -j 8 -o TRIMMED.$(inputs.readFilesIn1.basename)
  - position: 3
    shellQuote: false
    valueFrom: >-
      && cutadapt -j 8 -o TRIMMED.$(inputs.readFilesIn2.basename)


inputs:
  r1_max_len: {type: int, doc: "Set read1 to a fixed length. Useful for single cell experiments in which the number of cycles was improperly configured during sequencing. Recommend 28 for 10X v3 chemistry",
    inputBinding: { prefix: "--length", position: 1 } }
  r2_max_len: {type: int, doc: "Set read2 to a fixed length. Useful for single cell experiments in which the number of cycles was improperly configured during sequencing. Recommend 91 for 10X v3 chemistry",
    inputBinding: { prefix: "--length", position: 3 } }
  readFilesIn1: { type: File, doc: "read1 fastq file", inputBinding: {position: 2} }
  readFilesIn2: { type: File, doc: "read2 fastq file", inputBinding: {position: 4} }

outputs:
  trimmedReadsR1: { type: File, outputBinding: { glob: $("*TRIMMED." + inputs.readFilesIn1.basename) } }
  trimmedReadsR2: { type: File, outputBinding: { glob: $("*TRIMMED." + inputs.readFilesIn2.basename) } }
