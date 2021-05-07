cwlVersion: v1.0
class: CommandLineTool
id: merge_rnaseqc_results
doc: "Merge RNAseQC metrics files into one TSV file"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'ubuntu:20.04'
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: join_rec.sh
        entry:
          $include: ../scripts/join_rec.sh

baseCommand: []

arguments:
  - position: 1
    shellQuote: false
    valueFrom: |-
     qcs="${
       var arr = [];
       for (var i=0; i<inputs.metrics_files.length; i++)
        arr = arr.concat(inputs.metrics_files[i].path)
       return (arr.join(' '))
     }"
     /bin/bash join_rec.sh $qcs > $(inputs.output_basename).tsv

inputs:
  output_basename: {type: string, doc: "basename of output file"}
  metrics_files: {type: 'File[]', doc: "List of metrics TSVs to combine"}

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_basename).tsv
    doc: "TSV file with combined metrics"
