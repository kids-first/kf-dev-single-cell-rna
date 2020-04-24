cwlVersion: v1.0
class: CommandLineTool

id: cellranger_count

requirements:
  - class: DockerRequirement
    dockerPull: 'kfdrc/cellranger:3.1.0'
  - class: ResourceRequirement
    ramMin: 20000

baseCommand: [tar]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
     -xzf $(inputs.fastqs.path).tar.gz
     && tar -xzf $(inputs.fastqs.path).tar.gz
     && cellranger count --id=$(inputs.run_id) --fastqs=$(inputs.fastqs) --sample_name=$(inputs.sample_name) --reference=$(inputs.reference)
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
