cwlVersion: v1.2
class: CommandLineTool
id: hisat2_build_index
doc: "Build HiSAT2 Index. If gtf or snps added - will be an HGFM reference. Else if just fasta, will be a HFM reference"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/hisat2:2.2.1'
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
  ram: {type: 'int?', default: 200, doc: "In GB. Use 32 GB if just building a HFM reference"}
  cpus: {type: 'int?', default: 16, doc: "Number of CPUs to request",
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
  # advanced - if auto params cause memory errors
  noauto: { type: 'boolean?', default: false, doc: "Disable the default behavior whereby hisat2-build automatically selects values for the --bmax, --dcv and [--packed] parameters according to available memory",
    inputBinding: { prefix: "--noauto", position: 3 } }
  bmax: { type: 'int?', doc: "The maximum number of suffixes allowed in a block. Allowing more suffixes per block makes indexing faster, but increases peak memory usage. Setting this option overrides any previous setting for --bmax, or --bmaxdivn. Default (in terms of the --bmaxdivn parameter) is --bmaxdivn 4",
    inputBinding: { prefix: "--bmax", position: 3 } }
  bmaxdivn: { type: 'int?', doc: "The maximum number of suffixes allowed in a block, expressed as a fraction of the length of the reference. Setting this option overrides any previous setting for --bmax, or --bmaxdivn. Default: --bmaxdivn 4",
    inputBinding: { prefix: "--bmaxdivn", position: 3 } }
  dcv: { type: 'int?', doc: "Use <int> as the period for the difference-cover sample. A larger period yields less memory overhead, but may make suffix sorting slower, especially if repeats are present. Must be a power of 2 no greater than 4096. Default: 1024",
    inputBinding: { prefix: "--dcv", position: 3 } }
  nodc: { type: 'boolean?', doc: "Disable use of the difference-cover sample. Suffix sorting becomes quadratic-time in the worst case (where the worst case is an extremely repetitive reference). Default: off",
    inputBinding: { prefix: "--nodc", position: 3 } }

outputs:
  hisat_trans_tar:
    type: File
    outputBinding:
      glob: '*.tar.gz'
