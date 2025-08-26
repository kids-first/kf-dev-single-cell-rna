cwlVersion: v1.2
class: CommandLineTool
id: hisat2_format_gtf_ref
doc: "Convert gtf to HISAT2 format exon and splice files. Needed to build HGFM reference with transcripts"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/hisat2:2.2.1'
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 4
  - class: InlineJavascriptRequirement

baseCommand: []

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      hisat2_extract_splice_sites.py $(inputs.gtf.path) > $(inputs.output_basename).ss &&
      hisat2_extract_exons.py $(inputs.gtf.path) > $(inputs.output_basename).exon

inputs:
  gtf: {type: File, doc: "gtf reference file"}
  output_basename: {type: string, doc: "Output file basename. Recommend hisat2_gencodeNN"}

outputs:
  hisat2_ss:
    type: File
    outputBinding:
      glob: $(inputs.output_basename).ss
  hisat2_exon:
    type: File
    outputBinding:
      glob: $(inputs.output_basename).exon
