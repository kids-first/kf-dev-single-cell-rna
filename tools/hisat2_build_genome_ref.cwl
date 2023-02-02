cwlVersion: v1.2
class: CommandLineTool
id: hisat2_build_genome_ref
doc: "Build HiSAT2 Genome Reference"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0'
  - class: ResourceRequirement
    ramMin: ${return inputs.ram * 1000}
    coresMin: $(inputs.cpus)
  - class: InlineJavascriptRequirement

baseCommand: []
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      mkdir
  - position: 2
    shellQuote: false
    valueFrom: >-
      && hisat2-build
  - position: 4
    shellQuote: false
    valueFrom: >-
      $(inputs.output_basename)/$(inputs.output_basename) &&
      tar -czf $(inputs.output_basename).tar.gz $(inputs.output_basename)

inputs:
  fasta: { type: File, doc: "Reference fasta file",
    inputBinding: { position: 3 }}
  output_basename: {type: string, doc: "Output file basename. Recommend version_genome",
    inputBinding: { position: 1 }}
  ram: {type: ['null', int], default: 32, doc: "In GB"}
  cpus: {type: ['null', int], default: 16, doc: "Number of CPUs to request",
    inputBinding: { prefix: "-p", position: 2} }

outputs:
  hisat_genome_tar:
    type: File
    outputBinding:
      glob: '*.tar.gz'
