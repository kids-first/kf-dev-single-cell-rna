cwlVersion: v1.2
class: Workflow
id: kf_single_cell_10x_refinement
label: "KFDRC Single Cell RNA 10x Refinement Workflow"
doc: |
  # 10X Refinement Workflow

  The workflow script that runs the tools is `workflows/kf_single_cell_10x_refinement.cwl`

  [SoupX](https://github.com/constantAmateur/SoupX) is used for subtraction of the RNA background.
  [scDblFinder](https://github.com/plger/scDblFinder) is used to score and predict doublets.
  Decontaminated outputs are aggregated using the [Seurat](https://satijalab.org/seurat/) R package from the Satija lab at the New York Genome Center.
  Original workflow design heavily contributed to by Erin Reichenbee of DBHi.

  ## Software

  - soupX 4.1.0
  - scDblFinder 1.8.0
  - Seurat 4.0.4

  ## Inputs
  ### multi-step
   - `output_basename`: basename used to name output files
   - `sample_name`: used as prefix for finding fastqs to analyze, e.g. 1k_PBMCs_TotalSeq_B_3p_LT_antibody if the names of the underlying fastqs are of the form 1k_PBMCs_TotalSeq_B_3p_LT_antibody_S1_L001_I1_001.fastq.gz, one per input fastq in the same order
  ### soupX
   - `cellranger_matrix_raw`: Raw feature matrix file from Cellranger
   - `cellranger_matrix_filtered`: Filtered feature matrix file from Cellranger
   - `cellranger_cluster`: CSV containing cluster information from Cellranger
  ### scDblFinder
   - `align_qc_rds`: Align QC file frm D3b 10X alignment workflow
   - `seurat_raw_rds`: Seurat raw rds file from D3b 10X alignment workflow

  ### Outputs
  - `soupx_rplots`: PDF R plot made by soupX 
  - `scdblfinder_plot`: PDF cluster plots generated by scDblFinder
  - `scdblfinder_doublets`: TSV containing scoring matrix with doublets marked by scDblFinder
  - `scdblfinder_summary`: Summary stats of number and percent of doublets in library
  - `decontam_matrix`: RDS file containing merged count matrix from Seurat 
  - `decontam_object`: RDS file containing merged Seurat R object
requirements:
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
inputs:
  # multi-step
  output_basename: {type: string, doc: "basename used to name output files"}
  cellranger_matrix_raw: {type: 'File', doc: "raw_feature_bc_matrix file from cellranger count"}
  cellranger_matrix_filtered: {type: 'File', doc: "filtered_feature_bc_matrix file from cellranger count"}
  cellranger_cluster: {type: 'File', doc: "clusters.csv file from cellranger count"}
  align_qc_rds: {type: File, doc: "Align QC file frm D3b 10X alignment workflow"}
  seurat_raw_rds: {type: File, doc: "Seurat raw rds file from D3b 10X alignment workflow"}
  sample_name: {type: 'string', doc: "used as prefix for finding fastqs to analyze, e.g. 1k_PBMCs_TotalSeq_B_3p_LT_antibody if the
      names of the underlying fastqs are of the form 1k_PBMCs_TotalSeq_B_3p_LT_antibody_S1_L001_I1_001.fastq.gz, one per input fastq
      in the same order"}
outputs:
  soupx_rplots: {type: File, outputSource: rename_rplots/renamed_file}
  soupx_rds: {type: File, outputSource: rename_soupx_rds/renamed_file}
  scdblfinder_doublets: {type: File, outputSource: rename_doublets/renamed_file}
  scdblfinder_plot: {type: File, outputSource: rename_scdblfinder_plot/renamed_file}
  scdblfinder_summary: {type: File, outputSource: rename_scdblfinder_summary/renamed_file}
  seurat_filtered_rds: {type: File, outputSource: seurat_filter/filtered_rds}
steps:
  soupx:
    run: ../tools/soupx.cwl
    in:
      raw_matrix: cellranger_matrix_raw
      filtered_matrix: cellranger_matrix_filtered
      cluster_file: cellranger_cluster
      sample_name: sample_name
    out: [decontaminated_matrix_dir, decontaminated_matrix_rds, rplots]
  rename_rplots:
    run: ../tools/rename_file.cwl
    in:
      in_file: soupx/rplots
      out_filename:
        source: output_basename
        valueFrom: $(self).soupx.rplots.pdf
    out: [renamed_file]
  rename_soupx_rds:
    run: ../tools/rename_file.cwl
    in:
      in_file: soupx/decontaminated_matrix_rds
      out_filename:
        source: output_basename
        valueFrom: $(self).soupx.decontaminated_matrix.rds
    out: [renamed_file]
  scdblfinder:
    run: ../tools/scdblfinder.cwl
    in:
      seurat_raw_object: seurat_raw_rds
      sample_name: sample_name
    out: [result_dir]
  rename_doublets:
    run: ../tools/rename_file.cwl
    in:
      in_file:
        source: [scdblfinder/result_dir, sample_name]
        valueFrom: |
          ${for(var i in self[0].listing){
            var file = self[0].listing[i];
            if (file.basename == "doublets_" + self[1] + ".tsv"){
              return file;
              }
            }
          }
      out_filename:
        source: output_basename
        valueFrom: $(self).scDblFinder.doublets.tsv
    out: [renamed_file]
  rename_scdblfinder_plot:
    run: ../tools/rename_file.cwl
    in:
      in_file:
        source: [scdblfinder/result_dir, sample_name]
        valueFrom: |
          ${for(var i in self[0].listing){
            var file = self[0].listing[i];
            if (file.basename == self[1] + "-Doublets_prediction.pdf"){
              return file;}
              }
          }
      out_filename:
        source: output_basename
        valueFrom: $(self).scDblFinder.prediction.pdf
    out: [renamed_file]
  rename_scdblfinder_summary:
    run: ../tools/rename_file.cwl
    in:
      in_file:
        source: [scdblfinder/result_dir, sample_name]
        valueFrom: |
          ${for(var i in self[0].listing){
            var file = self[0].listing[i];
            if (file.basename == "doublets_metrics_" + self[1] + ".tsv"){
              return file;
              }
            }
          }
      out_filename:
        source: output_basename
        valueFrom: $(self).scDblFinder.summary_metrics.tsv
    out: [renamed_file]
  seurat_filter:
    run: ../tools/seurat_filter.cwl
    in:
      scdblfinder_tsv:
        source: [scdblfinder/result_dir, sample_name]
        valueFrom: |
          ${for(var i in self[0].listing){
            var file = self[0].listing[i];
            if (file.basename == "doublets_" + self[1] + ".tsv"){
              return file;
              }
            }
          }
      seurat_qc_rds: align_qc_rds
      soupx_rds: soupx/decontaminated_matrix_rds
      output_filename:
        source: output_basename
        valueFrom: $(self).seurat_qc.filtered.rds
    out: [filtered_rds]
sbg:license: Apache License 2.0
sbg:publisher: KFDRC
$namespaces:
  sbg: https://sevenbridges.com
hints:
- class: 'sbg:maxNumberOfParallelInstances'
  value: 2
