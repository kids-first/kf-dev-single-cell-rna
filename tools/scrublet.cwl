cwlVersion: v1.2
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
  input_matrix: {type: 'File?', loadListing: deep_listing, doc: "Cell ranger filtered h5 count matrix dir if preferred or matrix dir not available",
    inputBinding: { prefix: "--matrix", position: 1}}
  input_matrix_dir: {type: 'Directory?', loadListing: deep_listing, doc: "Cell ranger filtered count matrix dir if preferred or h5 not available",
    inputBinding: { prefix: "--matrix", position: 1}}
  output_basename: {type: string, doc: "Output files basename"}
  expected_doublet_rate: {type: 'float?', default: 0.06, doc: "expected doublet rate, usually specific to the method; default 0.06 for 10X", inputBinding: {prefix: -e, position: 2}}
  doublet_score_threshold: {type: 'float?', default: 0.25, doc: "doublet cut-off, cells with greater scores will be labelled as doublets; must be between 0 and 1", inputBinding: {prefix: -s, position: 2}}
  count_min: {type: 'int?', default: 2, doc: "minimum expression count to retain a gene", inputBinding: {prefix: -c, position: 2}}
  cell_min: {type: 'int?', default: 3, doc: "minimum number of cells a gene must be in to be retained", inputBinding: {prefix: -l, position: 2}}
  min_gene_variability_pctl: {type: 'int?', default: 85, doc: "Keep the most highly variable genes (in the top min_gene_variability_pctl percentile), as measured by the v-statistic", inputBinding: {prefix: -g, position: 2}}
  n_prin_comps: {type: 'int?', default: 30, doc: "Number of PCs to use for clustering", inputBinding: {prefix: -p, position: 2}}
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
