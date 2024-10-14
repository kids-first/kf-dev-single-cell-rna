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
      cutadapt -j 8 --length $(inputs.fixed_length[0]) -o TRIMMED.$(inputs.readFilesIn1.basename)
  - position: 3
    shellQuote: false
    valueFrom: >-
      && cutadapt -j 8 --length $(inputs.fixed_length[1]) -o TRIMMED.$(inputs.readFilesIn2.basename)


inputs:
  fixed_length: { type: 'int[]', doc: "Shorten reads to a fixed length. Useful for specialized protocols in which the instrument will read beyond useable sequence" }
  readFilesIn1: { type: File, doc: "read1 fastq file", inputBinding: {position: 2} }
  readFilesIn2: { type: File, doc: "read2 fastq file", inputBinding: {position: 4} }

outputs:
  trimmedReadsR1: { type: File, outputBinding: { glob: $("*TRIMMED." + inputs.readFilesIn1.basename) } }
  trimmedReadsR2: { type: File, outputBinding: { glob: $("*TRIMMED." + inputs.readFilesIn2.basename) } }
