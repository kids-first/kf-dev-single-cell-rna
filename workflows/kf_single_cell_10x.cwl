cwlVersion: v1.2
class: Workflow
id: kf_single_cell_10x

doc: |-
  # KFDRC single cell RNA 10x workflow
  Workflow for aligning, counting, and aggregation of single-cell RNA data generated by 10X.

  ![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

  The workflow runs [cellranger count](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count),
  on fastq files generated by the 10x single cell RNA workflow / methodology.
  Cell ranger count performs alignment, barcode counting, and filtering.
  [SoupX](https://github.com/constantAmateur/SoupX) is used for subtraction of the RNA background
  [Scrublet](https://github.com/swolock/scrublet) is used to scord and predict doublets
  Decontaminated outputs are aggregated using the [Seurat](https://satijalab.org/seurat/) R package from the Satija lab at the New York Genome Center.
  Additionally, [scanpy](https://scanpy.readthedocs.io/en/stable/) and [DropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html) are used for reading and writing intermediary files.

  ### Caveats:
  1. The reference file can be downloaded from 10x, however, the directory
     name contains several periods that must be changed before creating the
     reference tarball.
  1. The sample name is used as a prefix to find fastq files contained within
     the fastq tarball. cellranger expects reads to be in the format
     `[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz`. If the sample
     name is not found in the beginning of the fastq files or there are other
     characters before the _S1 part of the fastq name, the fastq will not be
     included in the analysis and if no files are found with the sample name,
     the workflow will fail.

requirements:
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}

inputs:
  output_basename: {type: string, doc: "basename used to name output files"}
  fastq_dirs: {type: 'Directory[]', doc: "directories of fastqs being run, one from each sample or well"}
  sample_name: {type: 'string[]', doc: "used as prefix for finding fastqs to analyze, e.g. 1k_PBMCs_TotalSeq_B_3p_LT_antibody if the names of the underlying fastqs are of the form 1k_PBMCs_TotalSeq_B_3p_LT_antibody_S1_L001_I1_001.fastq.gz, one per input fastq in the same order"}
  reference: {type: 'Directory', doc: "directory of reference files"}
  expected_doublet_rate: {type: 'float?', default: 0.06, doc: "expected doublet rate, usually specific to the method; default 0.06 for 10X"}
  doublet_score_threshold: {type: 'float?', default: 0.25, doc: "doublet cut-off, cells with greater scores will be labelled as doublets; must be between 0 and 1"}
  count_min: {type: 'int?', default: 2, doc: "minimum expression count to retain a gene"}
  cell_min: {type: 'int?', default: 3, doc: "minimum number of cells a gene must be in to be retained"}
  min_gene_variability_pctl: {type: 'int?', default: 85, doc: "Keep the most highly variable genes (in the top min_gene_variability_pctl percentile), as measured by the v-statistic"}
  n_prin_comps: {type: 'int?', default: 30, doc: "Number of PCs to use for clustering"}
  ram: {type: 'int?', default: 16, doc: "In GB"}
  cpus: {type: 'int?', default: 1, doc: "Number of CPUs to request"}

outputs:
  count_summary: {type: 'File[]', outputSource: count/output_summary}
  bam_out: {type: 'File[]', outputSource: count/bam}
  merged_decontam_matrix: {type: 'File', outputSource: merge/merged_matrix}
  merged_decontam_object: {type: 'File', outputSource: merge/merged_object}
  doublet_histogram: {type: 'File[]', outputSource: scrublet/score_histogram}

steps:

  count:
    run: ../tools/cellranger_count.cwl
    scatter: [fastqs, sample_name]
    scatterMethod: dotproduct
    in:
      run_id: output_basename
      fastqs: fastq_dirs
      sample_name: sample_name
      reference: reference
      return_h5:
        valueFrom: ${return Boolean(true)}
    out: [filtered_matrix_out, raw_matrix_out, bam, output_summary, molecule_info, whole_output_dir, cluster_file]

  soupx:
    run: ../tools/soupx.cwl
    scatter: [raw_matrix, filtered_matrix, sample_name, cluster_file]
    scatterMethod: dotproduct
    in:
      raw_matrix: count/raw_matrix_out
      filtered_matrix: count/filtered_matrix_out
      cluster_file: count/cluster_file
      sample_name: sample_name
    out: [decontaminated_matrix]

  scrublet:
    run: ../tools/scrublet.cwl
    scatter: [input_matrix, output_basename]
    scatterMethod: dotproduct
    in:
      input_matrix: soupx/decontaminated_matrix
      output_basename: sample_name
      expected_doublet_rate: expected_doublet_rate
      doublet_score_threshold: doublet_score_threshold
      count_min: count_min
      cell_min: cell_min
      min_gene_variability_pctl: min_gene_variability_pctl
      n_prin_comps: n_prin_comps
      ram: ram
      cpus: cpus
    out: [score_histogram, doublets_file]

  merge:
    run: ../tools/seurat_merge.cwl
    in:
      matrix_dirs: soupx/decontaminated_matrix
      doublets_files: scrublet/doublets_file
      output_name: output_basename
    out: [merged_matrix, merged_object]
