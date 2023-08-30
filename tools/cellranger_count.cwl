cwlVersion: v1.2
class: CommandLineTool
id: cellranger_count
doc: "Run cellranger count on a set of fastq files"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "pgc-images.sbgenomics.com/d3b-bixu/cellranger:6.1.2"
  - class: ResourceRequirement
    ramMin: $(inputs.cr_instance_ram * 1000)
    coresMin: $(inputs.localcores)
  - class: InlineJavascriptRequirement

baseCommand: [cellranger, count]

arguments:
  - position: 2
    shellQuote: false
    valueFrom: >-
     1>&2
     && mkdir $(inputs.sample_name)
     && mv $(inputs.run_id)/outs $(inputs.sample_name)/

inputs:
  localcores: { type: 'int?', doc: "Num cores to use", default: 16,
    inputBinding: {position: 1, prefix: "--localcores=", separate: false } }
  cr_instance_ram: { type: 'int?', doc: 'Ram in GB to make available to cell ranger count step', default: 64}
  run_id: { type: string, doc: "run id, used as basename for output", 
    inputBinding: { position: 1, prefix: "--id=", separate: false } }
  fastqs: { type: Directory, loadListing: deep_listing, doc: "directory of fastqs being run", 
    inputBinding: { position: 1, prefix: "--fastqs=", separate: false } }
  sample_name: { type: string, doc: "sample name, used as prefix for finding fastqs to analyze", 
    inputBinding: { position: 1, prefix: "--sample=", separate: false, shellQuote: false } }
  reference: { type: Directory, loadListing: deep_listing, doc: "directory of reference files", 
    inputBinding: { position: 1, prefix: "--transcriptome=", separate: false } }
  no_bam: { type: 'boolean?', doc: "Set to skip generating bam output. Good to keep bam for troubleshooting, but adds to computation time",
    inputBinding: { position: 1, prefix: "--no-bam"} }
  return_h5: { type: 'boolean?', doc: "TRUE: return h5 files or FALSE: return tarred matrix directories?" }
  include_introns: { type: 'boolean?', doc: "Include intronic reads in count",
    inputBinding: { position: 1, prefix: "--include-introns" } }
  chemistry: { type: ['null', {type: enum, name: chemistry, symbols: ["auto","threeprime","fiveprime","SC3Pv2","SC3Pv3","SC3Pv3LT","SC3Pv3HT","SC5P-PE","SC5P-R2","SC3Pv1","ARC-v1"]}],
    default: "auto", doc: "Chemistry used. auto is usually best. See README for exceptions", inputBinding: { position: 1, prefix: "--chemistry" } }

outputs:
  filtered_matrix_out:
    type: File
    outputBinding:
      glob: $(inputs.sample_name)/outs/filtered_feature_bc_matrix.*
    doc: "File containing the filtered feature matrix"
  raw_matrix_out:
    type: File
    outputBinding:
      glob: $(inputs.sample_name)/outs/raw_feature_bc_matrix.*
    doc: "File containing the raw feature matrix"
  bam:
    type: 'File?'
    outputBinding:
      glob: $(inputs.run_id)/outs/*.bam
    doc: "The bam file that was generated"
  whole_output_dir:
    type: Directory
    outputBinding:
      loadListing: deep_listing
      glob: $(inputs.sample_name)
  cluster_file:
    type: File
    outputBinding:
      glob: $(inputs.sample_name)/outs/analysis/clustering/graphclust/clusters.csv
    doc: "Cluster file, will be used by SoupX"
