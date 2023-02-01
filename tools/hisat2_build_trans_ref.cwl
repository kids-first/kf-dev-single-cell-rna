cwlVersion: v1.2
class: CommandLineTool
id: hisat2_build_trans_ref
doc: "Build HiSAT2 Transcript Reference"

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
      hisat2_extract_splice_sites.py $(inputs.gtf.path) > genome.ss &&
      hisat2_extract_exons.py $(inputs.gtf.path) > genome.exon &&
      mkdir $(inputs.output_basename) && 
      hisat2-build -p $(inputs.cpus) --exon genome.exon --ss genome.ss $(inputs.fasta.path) $(inputs.output_basename)/$(inputs.output_basename) &&
      tar -czf $(inputs.output_basename).tar.gz $(inputs.output_basename)

inputs:
  fasta: {type: File, doc: "tarball of reference files"}
  gtf: {type: File, doc: "gtf reference file"}
  output_basename: {type: string, doc: "Output file basename. Recommend gencode_vNN_trans"}
  ram: {type: ['null', int], default: 72, doc: "In GB"}
  cpus: {type: ['null', int], default: 36, doc: "Number of CPUs to request"}

outputs:
  hisat_trans_tar:
    type: File
    outputBinding:
      glob: '*.tar.gz'
