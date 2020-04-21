baseCommand: [cellranger, count]

requirements:
  - class: DockerRequirement
    dockerPull: 'kfdrc/cellranger:3.1.0'

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-

      --id=$(inputs.run_id)
      --fastqs=$(inputs.fastqs)
      --sample=$(inputs.sample_name)
      --transcriptome=$(inputs.reference)

      tar -czf
      run_count.tar.gz
      run_count_$(inputs.run_id)/

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
      glob: "run_count.tar.gz"
