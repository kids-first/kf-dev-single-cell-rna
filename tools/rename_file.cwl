cwlVersion: v1.2
class: CommandLineTool
id: rename_file
doc: >-
  Tool to rename files

requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: $(inputs.out_filename)
        entry: $(inputs.in_file)

baseCommand: [echo, done]

inputs:
  in_file:
    type: File
    doc: "Input file to rename"
  out_filename:
    type: string
    doc: "New file name for in_file"

outputs:
  renamed_file:
    type: File
    outputBinding:
      glob: $(inputs.out_filename)
