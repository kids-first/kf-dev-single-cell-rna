cwlVersion: v1.0
class: CommandLineTool
id: cellranger_count
doc: "Run cellranger count on a set of fastq files"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/cellranger:6.0'
  - class: ResourceRequirement
    ramMin: $(inputs.cr_instance_ram * 1000)
    coresMin: $(inputs.cores)
  - class: InlineJavascriptRequirement

baseCommand: [cellranger, count]

arguments:
  - position: 2
    shellQuote: false
    valueFrom: >-
     1>&2 && ${
       var sp = inputs.run_id + "/outs/" + inputs.run_id + "." + inputs.sample_name;
       var cmd = "mv " + inputs.run_id + "/outs/molecule_info.h5 " + sp + ".molecule_info.h5 && ";
       cmd += "mv " + inputs.run_id + "/outs/web_summary.html " + sp + ".web_summary.html && ";
       if ( inputs.bam_out != null){
         cmd += "mv " + inputs.run_id + "/outs/possorted_genome_bam.bam " + sp + ".bam && ";
       }
       if (inputs.return_h5){
         cmd += "mv " + inputs.run_id + "/outs/filtered_feature_bc_matrix.h5 " + sp + ".filtered_feature_bc_matrix.h5 && ";
         cmd += "mv " + inputs.run_id + "/outs/raw_feature_bc_matrix.h5 " + sp + ".raw_feature_bc_matrix.h5";
       }
       else {
         cmd += "mv " + inputs.run_id + "/outs/raw_feature_bc_matrix " + inputs.run_id + "/outs/raw_gene_bc_matrices && ";
         cmd += "mv " + inputs.run_id + "/outs/filtered_feature_bc_matrix " + inputs.run_id + "/outs/filtered_gene_bc_matrices && ";
         cmd += "mkdir " + inputs.run_id + "/outs/filtered_gene_bc_matrices/" + "GRCh38 && ";
         cmd += "mkdir " + inputs.run_id + "/outs/raw_gene_bc_matrices/" + "GRCh38 && ";
         cmd += "mv " + inputs.run_id + "/outs/filtered_gene_bc_matrices/*.gz " + inputs.run_id + "/outs/filtered_gene_bc_matrices/" + "GRCh38/ && ";
         cmd += "mv " + inputs.run_id + "/outs/raw_gene_bc_matrices/*.gz " + inputs.run_id + "/outs/raw_gene_bc_matrices/" + "GRCh38/ && ";
         cmd += "cp -r " + inputs.run_id + "/outs/filtered_gene_bc_matrices " + sp + ".filtered_feature_bc_matrix && ";
         cmd += "cp -r  " + inputs.run_id + "/outs/raw_gene_bc_matrices " + sp + ".raw_feature_bc_matrix && ";
         cmd += "tar -czf " + sp + ".filtered_feature_bc_matrix.tar.gz " + sp + ".filtered_feature_bc_matrix && ";
         cmd += "tar -czf " + sp + ".raw_feature_bc_matrix.tar.gz " + sp + ".raw_feature_bc_matrix";
       }
       return cmd;
     }

inputs:
  localcores: { type: 'int?', doc: "Num cores to use", default: 16,
    inputBinding: {position: 1, prefix: "--localcores", separate: false } }
  cr_instance_ram: { type: 'int?', doc: 'Ram in GB to make available to cell ranger count step', default: 64}
  run_id: { type: string, doc: "run id, used as basename for output", 
    inputBinding: { position: 1, prefix: "--id=", separate: false } }
  fastqs: { type: Directory, doc: "directory of fastqs being run", 
    inputBinding: { position: 1, prefix: "--fastqs=", separate: false } }
  sample_name: { type: string, doc: "sample name, used as prefix for finding fastqs to analyze", 
    inputBinding: { position: 1, prefix: "--sample=", separate: false } }
  reference: { type: Directory, doc: "directory of reference files", 
    inputBinding: { position: 1, prefix: "--transcriptome=", separate: false } }
  no_bam: { type: 'boolean?', doc: "Set to skip generating bam output. Good to keep bam for troubleshooting, but adds to computation time",
    inputBinding: { position: 1, prefix: "--no-bam"} }
  return_h5: { type: 'boolean?', doc: "TRUE: return h5 files or FALSE: return tarred matrix directories?" }

outputs:
  filtered_matrix_out:
    type: File
    outputBinding:
      glob: $(inputs.run_id)/outs/$(inputs.run_id).$(inputs.sample_name).filtered_feature_bc_matrix.*
    doc: "File containing the filtered feature matrix"
  raw_matrix_out:
    type: File
    outputBinding:
      glob: $(inputs.run_id)/outs/$(inputs.run_id).$(inputs.sample_name).raw_feature_bc_matrix.*
    doc: "File containing the raw feature matrix"
  bam:
    type: 'File?'
    outputBinding:
      glob: $(inputs.run_id)/outs/*.bam
    doc: "The bam file that was generated"
  output_summary:
    type: File
    outputBinding:
      glob: $(inputs.run_id)/outs/$(inputs.run_id).$(inputs.sample_name).web_summary.html
    doc: "HTML alignment summary"
  molecule_info:
    type: File
    outputBinding:
      glob: $(inputs.run_id)/outs/$(inputs.run_id).$(inputs.sample_name).molecule_info.h5
    doc: "Molecule info file, used by cellranger aggr"
  whole_output_dir:
    type: Directory
    outputBinding:
      glob: $(inputs.run_id)/outs
  cluster_file:
    type: File
    outputBinding:
      glob: $(inputs.run_id)/outs/analysis/clustering/graphclust/clusters.csv
    doc: "Cluster file, will be used by SoupX"
