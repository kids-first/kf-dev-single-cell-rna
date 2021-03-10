cwlVersion: v1.0
class: CommandLineTool
id: cellranger_aggr
doc: "Run cellranger aggr on a set of molecule_info files"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/cellranger:3.1'
  - class: ResourceRequirement
    ramMin: ${return inputs.ram * 1000}
    coresMin: $(inputs.cpus)
  - class: InlineJavascriptRequirement

baseCommand: [echo]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >
     "library_id,molecule_h5" > ./aggr_file.txt &&
     ${
       var cmd_str = "";
       for (var i=0; i<inputs.molecule_infos.length; i++){
         cmd_str += "echo " + inputs.molecule_infos[i].nameroot.concat(",", inputs.molecule_infos[i].path) + ">> ./aggr_file.txt;";
       }
       return (cmd_str);
     }

     cellranger aggr --localcores=$(inputs.cpus) --localmem=$(inputs.ram) --id=$(inputs.run_id) --csv=./aggr_file.txt &&

     ${
       var dir = inputs.run_id + "/outs/";
       var base_name = "./" + inputs.run_id + "_aggr_";
       var cmd = "mv " + dir + "web_summary.html " + base_name + "web_summary.html && ";
       if (inputs.return_h5){
         cmd += "mv " + dir + "filtered_feature_bc_matrix.h5 " + base_name + "filtered_feature_bc_matrix.h5 && ";
         cmd += "mv " + dir + "raw_feature_bc_matrix.h5 " + base_name + "raw_feature_bc_matrix.h5";
       }
       else {
         cmd += "mv " + dir + "filtered_feature_bc_matrix " + base_name + "filtered_feature_bc_matrix && ";
         cmd += "mv " + dir + "raw_feature_bc_matrix " + base_name + "raw_feature_bc_matrix && ";
         cmd += "tar -czf " + base_name + "filtered_feature_bc_matrix.tar.gz " + base_name + "filtered_feature_bc_matrix && ";
         cmd += "tar -czf " + base_name + "raw_feature_bc_matrix.tar.gz " + base_name + "raw_feature_bc_matrix";
       }
       return cmd;
     }

inputs:
  run_id: {type: string, doc: "run id, used as basename for output"}
  molecule_infos: {type: 'File[]', doc: "List of molecule info files to aggregate"}
  return_h5: {type: boolean?, doc: "Return h5 files or tarred matrix directories?"}
  ram: {type: ['null', int], default: 20, doc: "In GB"}
  cpus: {type: ['null', int], default: 16, doc: "Number of CPUs to request"}

outputs:
  filtered_matrix_out:
    type: File
    outputBinding:
      glob: ./$(inputs.run_id)_aggr_filtered_feature_bc_matrix.*
    doc: "Tarball or h5 file containing the filtered feature matrix"
  raw_matrix_out:
    type: File
    outputBinding:
      glob: ./$(inputs.run_id)_aggr_raw_feature_bc_matrix.*
    doc: "Tarball or h5 file containing the raw feature matrix"
  output_summary:
    type: File
    outputBinding:
      glob: ./$(inputs.run_id)_aggr_web_summary.html
    doc: "HTML alignment summary"
