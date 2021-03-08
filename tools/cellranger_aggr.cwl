cwlVersion: v1.0
class: CommandLineTool
id: cellranger_aggr
doc: "Run cellranger aggr on a set of molecule_info files"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/cellranger:3.1'
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 4
  - class: InlineJavascriptRequirement

baseCommand: [echo]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >
     "library_id,molecule_h5" > aggr_file.txt
     && aggr_data="${
       var arr = []
       for (var i=0; i<inputs.molecule_infos.length; i++)
        arr = arr.concat(inputs.molecule_infos[i].nameroot.concat(",", inputs.molecule_infos[i].path))
       return (arr.join('\n'))
     }"
     && echo -e $aggr_data >> aggr_file.txt

     cellranger aggr --localcores=4 --localmem=8 --id=$(inputs.run_id) --csv=aggr_file.txt &&

     ${
       var sp = inputs.run_id + "/outs/" + inputs.run_id + ".aggr";
       var cmd = "mv " + inputs.run_id + "/outs/web_summary.html " + sp + ".web_summary.html && ";
       if (inputs.return_h5){
         cmd += "mv " + inputs.run_id + "/outs/filtered_feature_bc_matrix.h5 " + sp + ".filtered_feature_bc_matrix.h5 && ";
         cmd += "mv " + inputs.run_id + "/outs/raw_feature_bc_matrix.h5 " + sp + ".raw_feature_bc_matrix.h5";
       }
       else {
         cmd += "mv " + inputs.run_id + "/outs/filtered_feature_bc_matrix " + sp + ".filtered_feature_bc_matrix &&" ;
         cmd += "mv " + inputs.run_id + "/outs/raw_feature_bc_matrix " + sp + ".raw_feature_bc_matrix && ";
         cmd += "tar -czf " + sp + ".filtered_feature_bc_matrix.tar.gz " + sp + ".filtered_feature_bc_matrix && ";
         cmd += "tar -czf " + sp + ".raw_feature_bc_matrix.tar.gz " + sp + ".raw_feature_bc_matrix";
       }
       return cmd;
     }

inputs:
  run_id: {type: string, doc: "run id, used as basename for output"}
  molecule_infos: {type: 'File[]', doc: "List of molecule info files to aggregate"}
  return_h5: {type: boolean?, doc: "Return h5 files or tarred matrix directories?"}

outputs:
  filtered_matrix_out:
    type: File
    outputBinding:
      glob: $(inputs.run_id)/outs/$(inputs.run_id).aggr.filtered_feature_bc_matrix.*
    doc: "Tarball or h5 file containing the filtered feature matrix"
  raw_matrix_out:
    type: File
    outputBinding:
      glob: $(inputs.run_id)/outs/$(inputs.run_id).aggr.raw_feature_bc_matrix.*
    doc: "Tarball or h5 file containing the raw feature matrix"
  output_summary:
    type: File
    outputBinding:
      glob: $(inputs.run_id)/outs/$(inputs.run_id).aggr.web_summary.html
    doc: "HTML alignment summary"
