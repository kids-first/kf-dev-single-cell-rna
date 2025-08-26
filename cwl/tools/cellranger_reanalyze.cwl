cwlVersion: v1.2
class: CommandLineTool
id: cellranger_reanalyze
doc: |-
  Run cellranger reanalyze
  [Reanalysis documentation from 10x](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/reanalyze)

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/cellranger:6.0'
  - class: ResourceRequirement
    ramMin: ${return inputs.ram * 1000}
    coresMin: $(inputs.cpus)
  - class: InlineJavascriptRequirement

baseCommand: [cellranger, reanalyze]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >
     &&
     ${
       var dir = inputs.run_id + "/outs/";
       var base_name = "./" + inputs.run_id + "_reanalysis_";
       var cmd = "mv " + dir + "web_summary.html " + base_name + "web_summary.html && ";
       if (inputs.return_h5){
         cmd += "mv " + dir + "filtered_feature_bc_matrix.h5 " + base_name + "filtered_feature_bc_matrix.h5";
       }
       else {
         cmd += "mv " + dir + "filtered_feature_bc_matrix " + base_name + "filtered_feature_bc_matrix && ";
         cmd += "tar -czf " + base_name + "filtered_feature_bc_matrix.tar.gz " + base_name + "filtered_feature_bc_matrix";
       }
       return cmd;
     }

inputs:
  run_id: {type: string, inputBinding: {prefix: --id}, doc: "run id, used as basename for output"}
  matrix: {type: 'File', inputBinding: {prefix: --matrix}, doc: "filtered or raw matrix file in h5 format"}
  extra_params: {type: 'File?', inputBinding: {prefix: --params}, doc: "optional csv file with a list of parameters to change"}
  agg_file: {type: 'File?', inputBinding: {prefix: --agg}, doc: "file used to aggregate datasets with cellranger aggr"}
  barcodes: {type: 'File?', inputBinding: {prefix: --barcodes}, doc: "csv file with a list of barcodes to use for analysis"}
  genes: {type: 'File?', inputBinding: {prefix: --genes}, doc: "file with gene ids to analyze"}
  exclude_genes: {type: 'File?', inputBinding: {prefix: --exclude-genes}, doc: "file with gene ids to exclude from analysis"}
  force_cells: {type: 'int?', inputBinding: {prefix: --force-cells}, doc: "force cellranger to use this number of cells instead of cellranger estimating them"}
  ram: {type: ['null', int], default: 20, inputBinding: {prefix: --localmem}, doc: "In GB"}
  cpus: {type: ['null', int], default: 16, inputBinding: {prefix: --localcores}, doc: "Number of CPUs to request"}
  return_h5: {type: 'boolean?', doc: "Return h5 files or tarred matrix directories?"}

outputs:
  filtered_matrix_out:
    type: File
    outputBinding:
      glob: $(inputs.run_id)_reanalysis_filtered_feature_bc_matrix.*
    doc: "Tarball or h5 file containing the filtered feature matrix"
  output_summary:
    type: File
    outputBinding:
      glob: $(inputs.run_id)_reanalysis_web_summary.html
    doc: "HTML alignment summary"
