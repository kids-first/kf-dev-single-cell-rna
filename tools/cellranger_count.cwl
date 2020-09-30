cwlVersion: v1.0
class: CommandLineTool
id: cellranger_count
doc: "Run cellranger count on a set of fastq files"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/cellranger:3.1.0'
  - class: ResourceRequirement
    ramMin: 20000
    coresMin: 16
  - class: InlineJavascriptRequirement

baseCommand: [tar, -xzf]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
     $(inputs.fastqs.path)
     && tar -xzf $(inputs.reference.path)
     && cellranger count --localcores=16 --id=$(inputs.run_id) --fastqs=./$(inputs.fastqs.nameroot.split('.')[0]) --sample=$(inputs.sample_name) --transcriptome=./$(inputs.reference.nameroot.split('.')[0])
     && mv $(inputs.run_id)/outs/molecule_info.h5 $(inputs.run_id)/outs/$(inputs.run_id).$(inputs.sample_name).molecule_info.h5
     && mv $(inputs.run_id)/outs/filtered_feature_bc_matrix $(inputs.run_id)/outs/$(inputs.run_id).$(inputs.sample_name).filtered_feature_bc_matrix
     && tar -C $(inputs.run_id)/outs/ -czf $(inputs.run_id).$(inputs.sample_name).filtered_feature_bc_matrix.tar.gz $(inputs.run_id).$(inputs.sample_name).filtered_feature_bc_matrix
     && mv $(inputs.run_id)/outs/raw_feature_bc_matrix $(inputs.run_id)/outs/$(inputs.run_id).$(inputs.sample_name).raw_feature_bc_matrix
     && tar -C $(inputs.run_id)/outs/ -czf $(inputs.run_id).$(inputs.sample_name).raw_feature_bc_matrix.tar.gz $(inputs.run_id).$(inputs.sample_name).raw_feature_bc_matrix

inputs:
  run_id: {type: string, doc: "run id, used as basename for output"}
  fastqs: {type: File, doc: "set of fastqs being run"}
  sample_name: {type: string, doc: "sample name, used as prefix for finding fastqs to analyze"}
  reference: {type: File, doc: "tarball of reference files"}
  return_h5: {type: ,doc: "Return h5 files or tarred matrix directories?"}

outputs:
  filtered_matrix_out:
    type: File
    outputBinding:
      glob: $(inputs.run_id).$(inputs.sample_name).filtered_feature_bc_matrix.tar.gz
    doc: "Tarball containing the filtered feature matrix"
  raw_matrix_out:
    type: File
    outputBinding:
      glob: $(inputs.run_id).$(inputs.sample_name).raw_feature_bc_matrix.tar.gz
    doc: "Tarball containing the raw feature matrix"
  bam:
    type: File
    outputBinding:
      glob: $(inputs.run_id)/outs/*.bam
    doc: "The bam file that was generated"
  output_summary:
    type: File
    outputBinding:
      glob: $(inputs.run_id)/outs/web_summary.html
    doc: "HTML alignment summary"
  molecule_info:
    type: File
    outputBinding:
      glob: $(inputs.run_id)/outs/$(inputs.run_id).$(inputs.sample_name).molecule_info.h5
    doc: "Molecule info file, used by cellranger aggr"
