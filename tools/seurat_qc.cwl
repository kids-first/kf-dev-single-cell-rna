cwlVersion: v1.2
class: CommandLineTool
id: seurat-qc
doc: "Run custom QC on 10X output"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "pgc-images.sbgenomics.com/brownm28/scrna_qc:v1.0.1"
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: seurat_alignment_qc.Rmd
        entry:
          $include: ../scripts/seurat_alignment_qc.Rmd
  - class: ResourceRequirement
    coresMin: 8
    ramMin: 16000

baseCommand: [Rscript, -e]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      "rmarkdown::render('seurat_alignment_qc.Rmd', clean = TRUE,
            params=list(scooter_path='/scooter', 
                        results_dir='.', 
                        data_path='$(inputs.filtered_bc_matrix_dir.path)', 
                        sample_name='$(inputs.sample_name)',
                        min_genes=$(inputs.min_genes),
                        max_genes=$(inputs.max_genes),
                        max_mt=$(inputs.max_mt),
                        normalize_method='$(inputs.normalize_method)',
                        num_pcs=$(inputs.num_pcs)
                        ))"
  - position: 2
    shellQuote: false
    valueFrom: >-
      && mv seurat_alignment_qc.nb.html Seurat_QC-$(inputs.sample_name).html
inputs:
  filtered_bc_matrix_dir: { type: Directory, loadListing: deep_listing }
  sample_name: { type: string }
  min_genes: { type: "int?", doc: "minimum number of genes per cell", default: 400 }
  max_genes: { type: "int?", doc: "maximum number of genes per cell", default: 4000 }
  max_mt: { type: "float?", doc: "maximum percent mitochondrial reads per cell", default: 5 }
  normalize_method: { type: ['null', {type: enum, name: normalize_method, symbols: ["log_norm","sct"]}],
    default: "log_norm", doc: "normalization method. One of log_norm or sct" }
  num_pcs: { type: "int?", doc: "number of PCs to calculate", default: 30 }

outputs:
  result_dir:
    type: Directory
    outputBinding:
      loadListing: deep_listing
      glob: 'Seurat_QC-$(inputs.sample_name)'
    doc: "Dir containing all files that were generated by Seurat"
  summary_html:
    type: File
    outputBinding:
      glob: "*.html"
    doc: "html file with qc results summary"
  rds:
    type: File
    outputBinding:
      glob: 'Seurat_QC-$(inputs.sample_name)/seurat_obj.rds'
    doc: "R data object for Seurat"
