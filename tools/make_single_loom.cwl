cwlVersion: v1.0
class: CommandLineTool
id: make_single_loom
doc: "Create a single sample loom file from the rsem gene output."

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/loompy:2.0.16'
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: make_rsem_loom.py
        entry:
          $include: ../scripts/make_rsem_loom.py

baseCommand: [python]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      make_rsem_loom.py --data $(inputs.rsem_matrix.path) --sample $(inputs.sample_name)

inputs:
  rsem_matrix: {type: 'File', doc: "Gene level output file from rsem"}
  sample_name: {type: string, doc: "Will be used as the sample's name inside the loom file and as the basename for the output loom file."}

outputs:
  loom_file:
    type: File
    outputBinding:
      glob: $(inputs.sample_name).loom
    doc: "Output loom file."
