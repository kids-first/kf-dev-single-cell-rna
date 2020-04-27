cwlVersion: v1.0
class: CommandLineTool

id: cellranger_count

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'kfdrc/cellranger:3.1.0'
  - class: ResourceRequirement
    ramMin: 20000
  - class: InlineJavascriptRequirement

baseCommand: [tar, -xzf]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
     $(inputs.fastqs.path)
     && tar -xzf $(inputs.reference.path)
     && cellranger count --id=$(inputs.run_id) --fastqs=./$(inputs.fastqs.nameroot.split('.')[0]) --sample=$(inputs.sample_name) --transcriptome=./$(inputs.reference.nameroot.split('.')[0])
     && tar -czf $(inputs.run_id).tar.gz $(inputs.run_id)

inputs:
  run_id:
    type: string
  fastqs:
    type: File
  sample_name:
    type: string
  reference:
    type: File

outputs:
  count_out:
    type: File
    outputBinding:
      glob: $(inputs.run_id).tar.gz
