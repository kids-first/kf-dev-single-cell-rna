cwlVersion: v1.2
class: Workflow
id: kf_single_cell_10x_refinement
label: "KFDRC Single Cell RNA 10x Refinement Workflow"
doc: |
  # 10X Refinement Workflow

  The workflow script that runs the tools is `workflows/kf_single_cell_10x_refinement.cwl`

  [SoupX](https://github.com/constantAmateur/SoupX) is used for subtraction of the RNA background.
  [Scrublet](https://github.com/swolock/scrublet) is used to score and predict doublets.
  Decontaminated outputs are aggregated using the [Seurat](https://satijalab.org/seurat/) R package from the Satija lab at the New York Genome Center.
  Original workflow design heavily contributed to by Erin Reichenbee of DBHi.

  ## Software

  - soupX 4.1.0
  - scrublet 0.2.3
  - Seurat 4.0.4

  ## Inputs
  ### multi-step
   - `output_basename`: basename used to name output files
   - `sample_name`: used as prefix for finding fastqs to analyze, e.g. 1k_PBMCs_TotalSeq_B_3p_LT_antibody if the names of the underlying fastqs are of the form 1k_PBMCs_TotalSeq_B_3p_LT_antibody_S1_L001_I1_001.fastq.gz, one per input fastq in the same order
  ### optional concat and rename step
   - `corrected_read_1_name`: corrected read one names in the 10x expected format 'SampleName_S1_L001_R1_001'. When provided, must be in the same order and same length as the sample name and corrected_read_2_name arrays.
   - `corrected_read_2_name`: corrected read two names in the 10x expected format 'SampleName_S1_L001_R2_001'. When provided, must be in the same order and same length as the sample name and corrected_read_1_name arrays.
  ### scrublet
   - `expected_doublet_rate`: expected doublet rate, usually specific to the method; default 0.06 for 10X
   - `doublet_score_threshold`: doublet cut-off, cells with greater scores will be labelled as doublets; must be between 0 and 1
   - `count_min`: minimum expression count to retain a gene
   - `cell_min`: minimum number of cells a gene must be in to be retained
   - `min_gene_variability_pctl`: Keep the most highly variable genes (in the top min_gene_variability_pctl percentile), as measured by the v-statistic
   - `n_prin_comps`: Number of PCs to use for clustering
   - `ram`: In GB
   - `cpus`: Number of CPUs to request

  ### Outputs
  - `soupx_rplots`: PDF R plot made by soupX
  - `scrublet_histogram`: PNG histogram made by scrublet
  - `scrublet_doublets`: CSV containing expression matrix with doublets removed by scrublet
  - `decontam_matrix`: RDS file containing merged count matrix from Seurat
  - `decontam_object`: RDS file containing merged Seurat R object
requirements:
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
inputs:
  # multi-step
  output_basename: {type: string, doc: "basename used to name output files"}
  cellranger_matrix_raw: {type: 'File', doc: "raw_feature_bc_matrix file from cellranger
      count"}
  cellranger_matrix_filtered: {type: 'File', doc: "filtered_feature_bc_matrix file
      from cellranger count"}
  cellranger_cluster: {type: 'File', doc: "clusters.csv file from cellranger count"}
  sample_name: {type: 'string', doc: "used as prefix for finding fastqs to analyze,
      e.g. 1k_PBMCs_TotalSeq_B_3p_LT_antibody if the names of the underlying fastqs
      are of the form 1k_PBMCs_TotalSeq_B_3p_LT_antibody_S1_L001_I1_001.fastq.gz,
      one per input fastq in the same order"}
  expected_doublet_rate: {type: 'float?', default: 0.06, doc: "expected doublet rate,
      usually specific to the method; default 0.06 for 10X"}
  doublet_score_threshold: {type: 'float?', default: 0.25, doc: "doublet cut-off,
      cells with greater scores will be labelled as doublets; must be between 0 and
      1"}
  count_min: {type: 'int?', default: 2, doc: "minimum expression count to retain a
      gene"}
  cell_min: {type: 'int?', default: 3, doc: "minimum number of cells a gene must be
      in to be retained"}
  min_gene_variability_pctl: {type: 'int?', default: 85, doc: "Keep the most highly
      variable genes (in the top min_gene_variability_pctl percentile), as measured
      by the v-statistic"}
  n_prin_comps: {type: 'int?', default: 10, doc: "Number of PCs to use for clustering"}
  scrublet_cpus: {type: 'int?', default: 1, doc: "Number of CPUs to request"}
  scrublet_ram: {type: 'int?', default: 16, doc: "In GB"}
outputs:
  soupx_rplots: {type: File, outputSource: rename_rplots/renamed_file}
  soupx_rds: { type: File, outputSource: rename_soupx_rds/renamed_file }
  scrublet_doublets: {type: File, outputSource: rename_doublets/renamed_file}
  scrublet_histogram: {type: File, outputSource: rename_histogram/renamed_file}
  decontam_matrix: {type: File, outputSource: seurat_merge/merged_matrix}
  decontam_object: {type: File, outputSource: seurat_merge/merged_object}
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
  scrublet:
    run: ../tools/scrublet.cwl
    in:
      input_matrix: cellranger_matrix_filtered
      output_basename: sample_name
      expected_doublet_rate: expected_doublet_rate
      doublet_score_threshold: doublet_score_threshold
      count_min: count_min
      cell_min: cell_min
      min_gene_variability_pctl: min_gene_variability_pctl
      n_prin_comps: n_prin_comps
      ram: scrublet_ram
      cpus: scrublet_cpus
    out: [score_histogram, doublets_file]
  rename_doublets:
    run: ../tools/rename_file.cwl
    in:
      in_file: scrublet/doublets_file
      out_filename:
        source: output_basename
        valueFrom: $(self).scrublet.doublets.csv
    out: [renamed_file]
  rename_histogram:
    run: ../tools/rename_file.cwl
    in:
      in_file: scrublet/score_histogram
      out_filename:
        source: output_basename
        valueFrom: $(self).scrublet.hist.png
    out: [renamed_file]
  seurat_merge:
    run: ../tools/seurat_merge.cwl
    in:
      matrix_dirs:
        source: soupx/decontaminated_matrix_dir
        valueFrom: |
          $([self])
      doublets_files:
        source: scrublet/doublets_file
        valueFrom: |
          $([self])
      output_name: output_basename
    out: [merged_matrix, merged_object]
sbg:license: Apache License 2.0
sbg:publisher: KFDRC
$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 2
