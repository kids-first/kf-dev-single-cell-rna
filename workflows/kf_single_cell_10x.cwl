cwlVersion: v1.0
class: Workflow
id: kf_singe_cell_10x
doc: |-
  # KFDRC single cell RNA 10x workflow
  Workflow for aligning and counting single cell RNA reads generated by 10x.

  ![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

  The workflow runs cellranger count and velocyto on fastq files generated by the
  10x single cell RNA workflow / methodology. Cell ranger count performs alignment,
  barcode counting, and filtering. Velocyto is used to compute RNA velocity using
  the counts of preprocessed and processed mRNA.

  ### Caveats:
  1. The fastqs and references inputs must be tarballs of folders containing the
     relevant files. The repeats and genes files are gtfs containg either repeats
     or genes and can be zipped.
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
  1. Multiple sample names can be given as a comma separated list.

inputs:
  run_id: {type: string, doc: "run id, used as basename cellranger for output"}
  fastqs: {type: File, doc: "tarball of fastqs being run"}
  sample_name: {type: string, doc: "sample name, used as prefix for finding fastqs to analyze"}
  reference: {type: File, doc: "tarball of reference files"}
  repeats: {type: File, doc: ".gtf file containing intervals to mask"}
  output_folder: {type: string, doc: "output folder"}
  genes: {type: File, doc: ".gtf file with genes to analyze"}

outputs:
  count_out: {type: File, outputSource: cellranger/count_out}
  velocyto_out: {type: File, outputSource: velocyto/velocyto_out}

steps:
  cellranger:
    run: ../tools/cellranger_count.cwl
    in:
      run_id: run_id
      fastqs: fastqs
      sample_name: sample_name
      reference: reference
    out: [count_out, bam, barcodes]

  velocyto:
    run: ../tools/velocyto.cwl
    in:
      sample_name: sample_name
      barcodes: cellranger/barcodes
      repeats: repeats
      output_folder: output_folder
      bam: cellranger/bam
      genes: genes
    out: [velocyto_out]
