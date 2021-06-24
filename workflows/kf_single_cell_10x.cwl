cwlVersion: v1.0
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
  count outputs will be aggregated using the [Seurat](https://satijalab.org/seurat/) R package from the Satija lab at the New York Genome Center.

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
  StepInputExpressionRequirement: {}

inputs:
  output_basename: {type: string, doc: "basename used to name output tarballs"}
  fastqs_tar: {type: 'File[]', doc: "tarball(s) of fastqs being run, one from each sample or well"}
  sample_name: {type: 'string[]', doc: "used as prefix for finding fastqs to analyze"}
  reference: {type: File, doc: "tarball of reference files"}
  count_h5_output: {type: boolean?, default: False, doc: "count returns hdf5 files [False]"}

outputs:
  count_summary: {type: 'File[]', outputSource: count/output_summary}
  bam_out: {type: 'File[]', outputSource: count/bam}
  decontam_matrix: {type: 'File[]', outputSource: soupx/decontaminated_matrix}

steps:

  count:
    run: ../tools/cellranger_count.cwl
    scatter: [fastqs, sample_name]
    scatterMethod: dotproduct
    in:
      run_id: output_basename
      fastqs: fastqs_tar
      sample_name: sample_name
      reference: reference
      return_h5: count_h5_output
    out: [filtered_matrix_out, raw_matrix_out, bam, output_summary, molecule_info, analysis]

  soupx:
    run: ../tools/soupx.cwl
    scatter: [count_dir, sample_name]
    scatterMethod: dotproduct
    in:
      count_dir:
        source: count/filtered_matrix_out
        valueFrom: |
          ${return {class: Directory, listing: [self.dirname]}}
      sample_name: sample_name
    out: [decontaminated_matrix]
