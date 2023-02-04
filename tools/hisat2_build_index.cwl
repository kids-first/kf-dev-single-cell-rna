cwlVersion: v1.2
class: CommandLineTool
id: hisat2_build_index
doc: "Build HiSAT2 Index. If gtf or snps added - will be an HGFM reference. Else if just fasta, will be a HFM reference"

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
      mkdir $(inputs.output_basename) &&
  - position: 2
    shellQuote: false
    valueFrom: >-
      hisat2-build
  - position: 99
    shellQuote: false
    valueFrom: >-
      $(inputs.output_basename)/$(inputs.output_basename) &&
      tar -czf $(inputs.output_basename).tar.gz $(inputs.output_basename)

inputs:
  ram: {type: ['null', int], default: 200, doc: "In GB. Use 32 GB if just building a HFM reference"}
  cpus: {type: ['null', int], default: 16, doc: "Number of CPUs to request",
    inputBinding: { prefix: "-p", position: 3 } }
  exon: {type: 'File?', doc: "HISAT2 format exon file, if adding transcripts",
    inputBinding: { prefix: "--exon", position: 4 } }
  splice: {type: 'File?', doc: "HISAT2 format splice (.ss) file, if adding transcripts",
    inputBinding: { prefix: "--ss", position: 4 } }
  snp: { type: 'File?', doc: "hisat2 format snp file if including snps",
    inputBinding: { prefix: "--snp", position: 5 } }
  haplotype: {type: 'File?', doc: "HISAT2 format haplotype (.haplotype) file, if including haplotypes",
    inputBinding: { prefix: "--haplotype", position: 5 } }
  large_index: { type: 'boolean?', doc: "If fasta >= 4BN bp or snp addition causes failure, make large index", default: false,
    inputBinding: { prefix: "--large-index", position: 6} }
  fasta: { type: File, doc: "reference fasta file",
    inputBinding: { position: 7 } }
  output_basename: {type: string, doc: "Output file basename. Recommend hisat2_gencodeNN_{snp and/or hap}_{trans/genome}"}

outputs:
  hisat_trans_tar:
    type: File
    outputBinding:
      glob: '*.tar.gz'
