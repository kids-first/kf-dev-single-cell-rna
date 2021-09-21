cwlVersion: v1.0
class: CommandLineTool
id: scrublet
doc: >-
  Run Scrublet on 10x output to remove predicted doublets
  https://github.com/swolock/scrublet

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/scrublet:0.2.3'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cpus)
    ramMin: ${return inputs.ram * 1000}
  - class: InitialWorkDirRequirement
    listing:
      - entryname: run_scrublet.py
        entry:
          $include: ../scripts/run_scrublet.py

baseCommand: [python]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      run_scrublet.py --matrix $(inputs.input_matrix.path) --output $(inputs.output_basename)

inputs:
  input_matrix: {type: Directory}
  output_basename: {type: string, doc: "Output files basename"}
  expected_doublet_rate: {type: 'float?', inputBinding: {prefix: -e}}
  doublet_score_threshold: {type: 'float?', inputBinding: {prefix: -s}}
  count_min: {type: 'int?', inputBinding: {prefix: -c}}
  cell_min: {type: 'int?', inputBinding: {prefix: -l}}
  min_gene_variability_pctl: {type: 'int?', inputBinding: {prefix: -g}}
  n_prin_comps: {type: 'int?', inputBinding: {prefix: -p}}
  ram: {type: 'int?', default: 16, doc: "In GB"}
  cpus: {type: 'int?', default: 1, doc: "Number of CPUs to request"}

outputs:
  score_histogram:
    type: File
    outputBinding:
      glob: $(inputs.output_basename).hist.png
    doc: "Histogram of doublet scores"
  doublets_file:
    type: File
    outputBinding:
      glob: $(inputs.output_basename).csv
    doc: "Expression matrix with doublets removed"
