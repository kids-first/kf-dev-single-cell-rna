cwlVersion: v1.0
class: Workflow
id: kf_singe_cell_10x
doc: |-
  # KFDRC single cell RNA 10x workflow
  Workflow for aligning and counting single cell RNA reads generated by 10x.

  ![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

  The workflow runs [cellranger count](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count),
  on fastq files generated by the 10x single cell RNA workflow / methodology.
  Cell ranger count performs alignment, barcode counting, and filtering.

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

inputs:
  output_basename: {type: string, doc: "basename used to name output tarballs"}
  fastqs: {type: File, doc: "tarball of fastqs being run"}
  sample_name: {type: string, doc: "used as prefix for finding fastqs to analyze"}
  reference: {type: File, doc: "tarball of reference files"}

outputs:
  summary_out: {type: File, outputSource: cellranger/output_summary}

steps:
  cellranger:
    run: ../tools/cellranger_count.cwl
    in:
      run_id: output_basename
      fastqs: fastqs
      sample_name: sample_name
      reference: reference
    out: [matrix_out, bam, barcodes, output_summary]
