cwlVersion: v1.2
class: CommandLineTool
id: soupx
doc: >-
  Run SoupX on 10x output to remove mRNA contamination
  https://github.com/constantAmateur/SoupX

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "pgc-images.sbgenomics.com/d3b-bixu/soupx_r:1.6.2"
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: run_soupX.R
        entry:
          $include: ../scripts/run_soupX.R

baseCommand: [Rscript]

arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      run_soupX.R

inputs:
  raw_matrix: { type: File, doc: "H5 file raw matrix produced by Cell Ranger count",
    inputBinding: { position: 1, prefix: "--raw" } }
  filtered_matrix: { type: File, doc: "H5 file filtered matrix produced by Cell Ranger count",
    inputBinding: { position: 1, prefix: "--filtered" } }
  sample_name: { type: string, doc: "Name of this sample",
    inputBinding: { position: 1, prefix: "--sample_name"} }
  results_dir: { type: 'string?',  doc: "Name of dir to create to store outputs", default: "./",
    inputBinding: { position: 1, prefix: "--results_dir" } }
  cluster_file: { type: 'File?',  doc: "CSV file clusters produced by Cell Ranger count. If other input or not available, leave out and tool will cluster",
    inputBinding: { position: 1, prefix: "--cluster" } }

outputs:
  decontaminated_matrix_rds:
    type: File
    outputBinding:
      glob: |
        $(inputs.results_dir == "./" ? inputs.sample_name + ".soupx.rds" : inputs.results_dir + "/" + inputs.sample_name + ".soupx.rds")
    doc: "rds containing decontaminated matrix"
  rplots:
    type: File
    outputBinding:
      glob: |
        $(inputs.results_dir == "./" ? inputs.sample_name + ".soupx.plots.pdf" : inputs.results_dir + "/" + inputs.sample_name + ".soupx.plots.pdf")
  marker_plots:
    type: File
    outputBinding:
      glob: |
        $(inputs.results_dir == "./" ? inputs.sample_name + ".soupx.plotMarkerDistribution.pdf" : inputs.results_dir + "/" + inputs.sample_name + ".soupx.plotMarkerDistribution.pdf")
