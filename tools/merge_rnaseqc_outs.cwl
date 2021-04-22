cwlVersion: v1.0
class: CommandLineTool
id: cellranger_aggr
doc: "Run cellranger aggr on a set of molecule_info files"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'ubuntu:20.04'
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: join_rec.sh
        entry: |-
          join_rec() {
              if [ $# -eq 1 ]; then
                join -1 1 -t $'\t' - "$1"
              else
                f=$1; shift
                join -1 1 -t $'\t' - "$f" | join_rec "$@"
              fi
          }
          if [ $# -le 2 ]; then
              join -1 1 -t $'\t' "$@"
          else
              f1=$1; f2=$2; shift 2
              join -1 1 -t $'\t' "$f1" "$f2" | join_rec "$@"
          fi

        writable: false

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
     echo $qcs &&
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
