cwlVersion: v1.0
class: CommandLineTool
id: soupx
doc: >-
  Run SoupX on 10x output to remove mRNA contamination
  https://github.com/constantAmateur/SoupX

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/soupx_r:4.1.0_SoupX'
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: run_soupX.R
        entry:
          $include: ../scripts/run_soupX.R

baseCommand: [Rscript]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      run_soupX.R --countmatrixdir $(inputs.count_dir.path) --sample_name $(inputs.sample_name)

inputs:
  count_dir:
    type: Directory
    doc: "Directory containing raw and filtered mtx files produced by Cell Ranger count"
  sample_name:
    type: string
    doc: "Name of this sample"

outputs:
  decontaminated_matrix:
    type: File
    outputBinding:
      glob: $(inputs.sample_name).decontam.RDS
    doc: "RDS file containing decontaminated count matrix"
