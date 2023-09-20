cwlVersion: v1.2
class: Workflow
id: kf_single_cell_10x_alignment
label: "KFDRC Single Cell RNA 10x Alignment Workflow"
doc: |
  # 10X Alignment Workflow

  The workflow script that runs the tools is `workflows/kf_single_cell_10x_alignment.cwl`

  The workflow runs [cellranger count](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/6.0/using/count),
  on fastq files generated by the 10x single cell RNA workflow methodology.
  Cell ranger count performs alignment, barcode counting, and filtering.
  A custom QC R markdown notebook developed by @AntoniaChroni is also run, which includes as it's main engine Seurat and [scooter](https://github.com/igordot/scooter)

  ## Software

  - Cellranger 6.1.2
  - Seurat 4.0.4

  ## Inputs
  ### multi-step
   - `output_basename`: basename used to name output files
   - `sample_name`: used as prefix for finding fastqs to analyze, e.g. 1k_PBMCs_TotalSeq_B_3p_LT_antibody if the names of the underlying fastqs are of the form 1k_PBMCs_TotalSeq_B_3p_LT_antibody_S1_L001_I1_001.fastq.gz, one per input fastq in the same order
  ### optional concat and rename step
   - `corrected_read_1_name`: corrected read one names in the 10x expected format 'SampleName_S1_L001_R1_001'. When provided, must be in the same order and same length as the sample name and corrected_read_2_name arrays.
   - `corrected_read_2_name`: corrected read two names in the 10x expected format 'SampleName_S1_L001_R2_001'. When provided, must be in the same order and same length as the sample name and corrected_read_1_name arrays.
  ### cell ranger
   - `cr_localcores`: Num cores to use for cell ranger, default: 36
   - `cr_instance_ram`: Ram in GB to make available to cell ranger count step, default: 64
   - `fastq_dir`: directory of fastqs being run. If formatting needed, use r1 and r2 fastqs input instead
   - `r1_fastqs`: If fastqs need to be concat from an old format, populate this
   - `r2_fastqs`: If fastqs need to be concat from an old format, populate this
   - `reference`: directory of reference files
   - `no_bam`: Set to skip generating bam output. Good to keep bam for troubleshooting, but adds to computation time
   - `chemistry`:
     - `auto`: for auto-detection (default)
     - `threeprime`: for Single Cell 3′
     - `fiveprime`: for Single Cell 5′
     - `SC3Pv2`: for Single Cell 3′ v2
     - `SC3Pv3`: for Single Cell 3′ v3
     - `SC3Pv3LT`: for Single Cell 3′ v3 LT
     - `SC3Pv3HT`: for Single Cell 3′ v3 HT
     - `SC5P-PE`: for Single Cell 5′ paired-end (both R1 and R2 are used for alignment)
     - `SC5P-R2`: for Single Cell 5′ R2-only (where only R2 is used for alignment)
     - `SC3Pv1`: for Single Cell 3′ v1. NOTE: this mode cannot be auto-detected. It must be set explicitly with this option
     - `ARC-v1`: for analyzing the GEX portion of multiome data. NOTE: this mode cannot be auto-detected
  ### seurat qc
   - `seurat_qc_min_genes`: minimum number of genes per cell
   - `seurat_qc_max_genes`: maximum number of genes per cell
   - `seurat_qc_max_mt`: maximum percent mitochondrial reads per cell
   - `seurat_qc_normalize_method`: normalization method. One of log_norm or sct
   - `seurat_qc_nfeatures`: number of variable features to extract
   - `seurat_qc_num_pcs`: number of PCs to calculate

  ### Outputs
  - `bam_out`: BAM generated by Cellranger Count
  - `debug_cr_file_outputs`: TAR.GZ file of the output directory produced by Cellranger Count
  - `seurat_qc_html`: HTML of QC metrics generated by Seurat
  - `seurat_qc_rds`: RDS object file generated by Seurat
requirements:
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
inputs:
  # multi-step
  output_basename: {type: string, doc: "basename used to name output files"}
  sample_name: {type: 'string', doc: "used as prefix for finding fastqs to analyze,
      e.g. 1k_PBMCs_TotalSeq_B_3p_LT_antibody if the names of the underlying fastqs
      are of the form 1k_PBMCs_TotalSeq_B_3p_LT_antibody_S1_L001_I1_001.fastq.gz,
      one per input fastq in the same order"}
  corrected_read_1_name: {type: 'string?', doc: "corrected read one names in the 10x
      expected format 'SampleName_S1_L001_R1_001'. When provided, must be in the same
      order and same length as the sample name and corrected_read_2_name arrays."}
  corrected_read_2_name: {type: 'string?', doc: "corrected read two names in the 10x
      expected format 'SampleName_S1_L001_R2_001'. When provided, must be in the same
      order and same length as the sample name and corrected_read_1_name arrays."}
  cr_localcores: {type: 'int?', doc: "Num cores to use for cell ranger", default: 36}
  cr_instance_ram: {type: 'int?', doc: 'Ram in GB to make available to cell ranger
      count step', default: 64}
  fastq_dir: {type: 'Directory?', loadListing: deep_listing, doc: "directory of fastqs
      being run. If formatting needed, use r1 and r2 fastqs input instead"}
  r1_fastqs: {type: 'File[]?', doc: "If fastqs need to be concat from an old format,
      populate this"}
  r2_fastqs: {type: 'File[]?', doc: "If fastqs need to be concat from an old format,
      populate this"}
  reference: {type: 'Directory', loadListing: deep_listing, doc: "directory of reference
      files"}
  no_bam: {type: 'boolean?', doc: "Set to skip generating bam output. Good to keep
      bam for troubleshooting, but adds to computation time"}
  include_introns: {type: 'boolean?', doc: "Include intronic reads in count", default: false}
  chemistry: {type: ['null', {type: enum, name: chemistry, symbols: ["auto", "threeprime",
          "fiveprime", "SC3Pv2", "SC3Pv3", "SC3Pv3LT", "SC3Pv3HT", "SC5P-PE", "SC5P-R2",
          "SC3Pv1", "ARC-v1"]}], default: "auto", doc: "Chemistry used. auto is usually
      best. See README for exceptions"}
  seurat_qc_min_genes: {type: "int?", doc: "minimum number of genes per cell", default: 400}
  seurat_qc_max_genes: {type: "int?", doc: "maximum number of genes per cell", default: 4000}
  seurat_qc_max_mt: {type: "int?", doc: "maximum percent mitochondrial reads per cell",
    default: 5}
  seurat_qc_normalize_method: {type: ['null', {type: enum, name: normalize_method,
        symbols: ["log_norm", "sct"]}], default: "log_norm", doc: "normalization method.
      One of log_norm or sct"}
  seurat_qc_nfeatures: {type: "int?", doc: "number of variable features to extract",
    default: 2000}
  seurat_qc_num_pcs: {type: "int?", doc: "number of PCs to calculate", default: 30}
outputs:
  debug_cr_file_outputs: {type: 'File', outputSource: tar_count_outdir/output }
  cellranger_bam: {type: 'File?', outputSource: rename_bam/renamed_file }
  cellranger_matrix_filtered: {type: 'File', outputSource: rename_matrix_filtered/renamed_file }
  cellranger_matrix_raw: {type: 'File', outputSource: rename_matrix_raw/renamed_file }
  cellranger_cluster: {type: 'File', outputSource: rename_clusters/renamed_file }
  seurat_qc_html: {type: File, outputSource: rename_seurat_html/renamed_file }
  seurat_qc_rds: {type: File, outputSource: rename_seurat_rds/renamed_file }
steps:
  concat_rename_fastq:
    run: ../tools/concat_rename_fastq.cwl
    when: $(inputs.r1_fastqs != null)
    in:
      r1_fastqs: r1_fastqs
      r2_fastqs: r2_fastqs
      sample_name: sample_name
      corrected_read_1_name: corrected_read_1_name
      corrected_read_2_name: corrected_read_2_name
    out: [renamed_dir]
  cellranger_count:
    run: ../tools/cellranger_count.cwl
    in:
      localcores: cr_localcores
      cr_instance_ram: cr_instance_ram
      run_id: output_basename
      fastqs:
        source: [concat_rename_fastq/renamed_dir, fastq_dir]
        pickValue: first_non_null
      sample_name: sample_name
      reference: reference
      no_bam: no_bam
      return_h5:
        valueFrom: ${return Boolean(true)}
      include_introns: include_introns
      chemistry: chemistry
    out: [filtered_matrix_out, raw_matrix_out, bam, whole_output_dir, cluster_file]
  rename_bam:
    run: ../tools/rename_file.cwl
    in:
      in_file: cellranger_count/bam
      out_filename:
        source: output_basename
        valueFrom: $(self).cellranger.count.possorted.genome.bam
    out: [renamed_file]
  rename_clusters:
    run: ../tools/rename_file.cwl
    in:
      in_file: cellranger_count/cluster_file
      out_filename:
        source: output_basename
        valueFrom: $(self).cellranger.count.clusters.csv
    out: [renamed_file]
  rename_matrix_filtered:
    run: ../tools/rename_file.cwl
    in:
      in_file: cellranger_count/filtered_matrix_out
      out_filename:
        source: output_basename
        valueFrom: $(self).cellranger.count.filtered_feature_bc_matrix.h5
    out: [renamed_file]
  rename_matrix_raw:
    run: ../tools/rename_file.cwl
    in:
      in_file: cellranger_count/raw_matrix_out
      out_filename:
        source: output_basename
        valueFrom: $(self).cellranger.count.raw_feature_bc_matrix.h5
    out: [renamed_file]
  tar_count_outdir:
    run: ../tools/tar.cwl
    in:
      output_filename:
        source: output_basename
        valueFrom: $(self).cellranger.count.tar.gz
      input_dir: cellranger_count/whole_output_dir
    out: [output]
  seurat_qc:
    run: ../tools/seurat_qc.cwl
    in:
      filtered_bc_matrix_dir: cellranger_count/whole_output_dir
      sample_name: sample_name
      min_genes: seurat_qc_min_genes
      max_genes: seurat_qc_max_genes
      max_mt: seurat_qc_max_mt
      normalize_method: seurat_qc_normalize_method
      nfeatures: seurat_qc_nfeatures
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
sbg:license: Apache License 2.0
sbg:publisher: KFDRC
$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 2
