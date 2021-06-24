cwlVersion: v1.2
class: CommandLineTool
id: loom_to_mtx
doc: "Convert a loom file to a Matrix Market file that Seurat can read."

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/loompy:2.0.16'
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: loom_to_mtx.py
        entry:
          $include: ../scripts/loom_to_mtx.py

baseCommand: [mkdir, -p]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.output_basename) &&
      python loom_to_mtx.py $(inputs.loom_file.path) $(inputs.output_basename) &&
      gzip $(inputs.output_basename)/* &&
      tar cvzf $(inputs.output_basename).tar.gz $(inputs.output_basename)

inputs:
  loom_file: {type: 'File', doc: "Loom file to convert"}
  output_basename: {type: string, doc: "Basename for the output tarball"}

outputs:
  output_tar:
    type: File
    outputBinding:
      glob: $(inputs.output_basename).tar.gz
    doc: "Output tarball"
