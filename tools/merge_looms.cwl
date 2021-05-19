cwlVersion: v1.0
class: CommandLineTool
id: merge_looms
doc: "Merge several loom files into one"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/loompy:2.0.16'
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: assemble_ss2_output.py
        entry:
          $include: ../scripts/assemble_ss2_output.py

baseCommand: []

arguments:
  - position: 1
    shellQuote: false
    valueFrom: |-
     locs="${
       var arr = [];
       for (var i=0; i<inputs.loom_files.length; i++)
        arr = arr.concat(inputs.loom_files[i].path)
       return (arr.join(' '))
     }"
     looms="${
       var arr = [];
       for (var i=0; i<inputs.loom_files.length; i++)
        arr = arr.concat("/tmp/" + inputs.loom_files[i].basename)
       return (arr.join(','))
     }"
     cp -t /tmp/ $locs
     python assemble_ss2_output.py --base $(inputs.output_basename) --files $looms

inputs:
  output_basename: {type: string, doc: "basename of output file"}
  loom_files: {type: 'File[]', doc: "List of loom files to combine"}

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_basename).loom
    doc: "Combined matrix loom file"
