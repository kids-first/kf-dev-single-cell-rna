cwlVersion: v1.0
class: CommandLineTool

id: cellranger_count

baseCommand: [cellranger, count]

requirements:
  - class: DockerRequirement
    dockerPull: 'kfdrc/cellranger:3.1.0'

inputs:
  run_id:
    type: string
    inputBinding:
      position: 1
      prefix: --id=
  fastqs:
    type: string
    inputBinding:
      position: 2
      prefix: --fastqs=
  sample_name:
    type: string
    inputBinding:
      position: 3
      prefix: --sample=
  reference:
    type: string
    inputBinding:
      position: 4
      prefix: --transcriptome=

outputs:
  count_out:
    type: Directory
    outputBinding:
      glob: "run_count_$(inputs.run_id)/outs"
