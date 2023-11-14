cwlVersion: v1.2
class: Workflow
id: kf-qc-repeat
label: Kids First DRC 10X QC Repeat Workflow
doc: "Repeat QC step using already-aligned data with new params"

requirements:
  - class: MultipleInputFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

inputs:
  debug_cr_file_outputs: {type: 'File', doc: "Tar gzipped outputs from cell ranger run" }
  output_basename: {type: string, doc: "basename used to name output files"}
  sample_name: {type: 'string', doc: "used as prefix for finding fastqs to analyze,
      e.g. 1k_PBMCs_TotalSeq_B_3p_LT_antibody if the names of the underlying fastqs
      are of the form 1k_PBMCs_TotalSeq_B_3p_LT_antibody_S1_L001_I1_001.fastq.gz,
      one per input fastq in the same order"}
  seurat_qc_min_genes: {type: "int?", doc: "minimum number of genes per cell", default: 400}
  seurat_qc_max_genes: {type: "int?", doc: "maximum number of genes per cell", default: 4000}
  seurat_qc_max_mt: {type: "float?", doc: "maximum percent mitochondrial reads per cell",
    default: 5}
  seurat_qc_normalize_method: {type: ['null', {type: enum, name: normalize_method,
        symbols: ["log_norm", "sct"]}], default: "log_norm", doc: "normalization method.
      One of log_norm or sct"}
  seurat_qc_num_pcs: {type: "int?", doc: "number of PCs to calculate", default: 30}

outputs:
  seurat_qc_html: {type: File, outputSource: rename_seurat_html/renamed_file }
  seurat_qc_rds: {type: File, outputSource: rename_seurat_rds/renamed_file }

steps:
  untar_dir:
    run: ../tools/untar_dir.cwl
    in:
      tarfile: debug_cr_file_outputs
      untar_dir_name: sample_name
    out: [outdir]
  seurat_qc:
    run: ../tools/seurat_qc.cwl
    in:
      filtered_bc_matrix_dir: untar_dir/outdir
      sample_name: sample_name
      min_genes: seurat_qc_min_genes
      max_genes: seurat_qc_max_genes
      max_mt: seurat_qc_max_mt
      normalize_method: seurat_qc_normalize_method
      num_pcs: seurat_qc_num_pcs
    out: [result_dir, summary_html, rds]
  rename_seurat_html:
    run: ../tools/rename_file.cwl
    in:
      in_file: seurat_qc/summary_html
      out_filename:
        source: output_basename
        valueFrom: $(self).seurat.qc.html
    out: [renamed_file]
  rename_seurat_rds:
    run: ../tools/rename_file.cwl
    in:
      in_file: seurat_qc/rds
      out_filename:
        source: output_basename
        valueFrom: $(self).seurat.qc.rds
    out: [renamed_file]