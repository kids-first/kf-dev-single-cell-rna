cwlVersion: v1.0
class: CommandLineTool

id: cellranger_count

baseCommand: ['cellranger', 'count']

requirements:
  - class: DockerRequirement
    dockerPull: 'kfdrc/cellranger:3.1.0'
  - class: ResourceRequirement
    ramMin: 8000

inputs:
  run_id:
    type: string
    inputBinding:
        position: 1
        prefix: --id=
        separate: false
  fastqs:
    type: Directory
    inputBinding:
        position: 2
        prefix: --fastqs=
        separate: false
  sample_name:
    type: string
    inputBinding:
        position: 3
        prefix: --sample=
        separate: false
  reference:
    type: Directory
    inputBinding:
        position: 4
        prefix: --transcriptome=
        separate: false

outputs:
  count_out:
    type: Directory
    outputBinding:
      glob: "$(inputs.run_id)"
