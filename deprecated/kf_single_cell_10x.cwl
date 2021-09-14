cwlVersion: v1.0
class: Workflow
id: kf_single_cell_10x

doc: |-
  # KFDRC single cell RNA 10x workflow
  Workflow for aligning, counting, aggregation and basic analysis of  single-cell RNA data generated by 10X.

  ![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

  The workflow runs [cellranger count](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count),
  on fastq files generated by the 10x single cell RNA workflow / methodology.
  Cell ranger count performs alignment, barcode counting, and filtering.
  count outputs are aggregated using [cellranger aggr](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate)
  The aggregated data is filtered and analyzed using the [Seurat](https://satijalab.org/seurat/) R package from the Satija lab at the New York Genome Center.
  Basic functionality includes removal of abnormal cells, dimensionality reduction, identification of differentially expressed features, and clustering.

  ### Caveats:
  1. The fastqs and references inputs must be tarballs of folders containing the
     relevant files.
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

inputs:
  output_basename: {type: string, doc: "basename used to name output tarballs"}
  fastqs_tar: {type: 'File[]', doc: "tarball(s) of fastqs being run, one from each sample or well"}
  sample_name: {type: string, doc: "used as prefix for finding fastqs to analyze"}
  reference: {type: File, doc: "tarball of reference files"}
  count_h5_output: {type: boolean?, default: False, doc: "count returns hdf5 files [False]"}
  aggr_h5_output: {type: boolean?, default: False, doc: "aggr returns hdf5 files [False]"}
  seurat_min_features: {type: int?, default: 200, doc: "Minimum number of genes observed in a cell to retain"}
  seurat_max_features: {type: int?, default: 2500, doc: "Maximum number of genes observed in a cell to retain"}
  seurat_max_mt: {type: int?, default: 5, doc: "Maximum mitochondrial percentage observed in a cell to retain"}
  seurat_norm_method: {type: string?, default: "LogNormalize", doc: "Normalization to apply to counts (LogNormalize, CLR, RC)"}
  seurat_retain_features: {type: int?, default: 2000, doc: "Number of most-variable features to initially retain"}
  seurat_nheatmap: {type: int?, default: 10, doc: "Number of principal components for which to produce heatmaps"}
  seurat_num_pcs: {type: int?, default: 10, doc: "Number of principal components to retain for clustering"}
  seurat_knn_granularity: {type: float?, default: 0.5, doc: "KNN clustering granularity parameter"}

outputs:
  count_summary: {type: 'File[]', outputSource: count/output_summary}
  aggr_summary: {type: File, outputSource: aggr/output_summary}
  filt_matrix_out: {type: File, outputSource: aggr/filtered_matrix_out}
  raw_matrix_out: {type: File, outputSource: aggr/raw_matrix_out}
  seurat_tarball: {type: File, outputSource: seurat/tarball}

steps:

  count:
    run: ../tools/cellranger_count.cwl
    scatter: [fastqs]
    in:
      run_id: output_basename
      fastqs: fastqs_tar
      sample_name: sample_name
      reference: reference
      return_h5: count_h5_output
    out: [filtered_matrix_out, raw_matrix_out, bam, output_summary, molecule_info]

  aggr:
    run: ../tools/cellranger_aggr.cwl
    in:
      run_id: output_basename
      molecule_infos: count/molecule_info
      return_h5: aggr_h5_output
    out: [filtered_matrix_out, raw_matrix_out, output_summary]

  seurat:
    run: ../tools/seurat.cwl
    in:
      name: sample_name
      scRNA_cts_tar: aggr/filtered_matrix_out
      output_basename: output_basename 
      min_features: seurat_min_features 
      max_features: seurat_max_features
      max_mt: seurat_max_mt 
      norm_method: seurat_norm_method
      retain_features: seurat_retain_features
      nheatmap: seurat_nheatmap
      num_pcs: seurat_num_pcs
      knn_granularity: seurat_knn_granularity
    out: [tarball]