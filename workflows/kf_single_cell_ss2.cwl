cwlVersion: v1.2
class: Workflow
id: kf_singe_cell_ss2
label: Kids First DRC single cell RNA Smart Seq 2 Workflow
doc: |-
  # KFDRC single cell RNA Smart Seq 2 workflow
  Workflow for aligning and counting single cell RNA reads generated by Smart Seq 2.

  ![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)
  The workflow uses [HISAT2](http://daehwankimlab.github.io/hisat2/) for alignment, [RNAseQC](https://github.com/getzlab/rnaseqc) to collect sequencing metrics, and [RSEM](https://deweylab.github.io/RSEM/) to calculate gene expression.
  The outputs of the workflow are a matrix of cells gene counts in a loom file, a tarball containing a matrix of the gene counts in Matrix Market format, and a tsv with collected sequencing metrics from all cells analyzed.
  Basic functionality includes removal of abnormal cells, dimensionality reduction, identification of differentially expressed features, and clustering.

  Data are aligned to both the genome and the transcriptome. Genome aligned data are used to collect sequencing metrics and transcriptome aligned data are used to calculate expression.

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

inputs:
  final_output_basename: {type: string, doc: "Output basename for workflow output files"}
  input_dir: {type: 'Directory', doc: "Directory containing fastq files"}
  hisat_genome_ref: {type: 'File', doc: "Hisat 2 genome reference"}
  hisat_trans_ref: {type: 'File', doc: "Hisat 2 transcriptome reference"}
  rnaseqc_gtf: {type: "File", doc: "gtf file used by RNAseQC", "sbg:suggestedValue": {class: 'File', path: '62853e7ad63f7c6d8d7ae5a3', name: 'gencode.v39.primary_assembly.rnaseqc.stranded.gtf'}}
  rsem_reference: {type: "File", doc: "RSEM reference file", "sbg:suggestedValue": {class: 'File', path: '62853e7ad63f7c6d8d7ae5a5', name: 'RSEM_GENCODE39.tar.gz'}}
  hisat_cpus: { type: 'int?', default: 4, doc: "CPUs to allocate to call task"}
  hisat_ram: { type: 'int?', default: 8, doc: "RAM to allocate to call task in gb"}
  rsem_cpus: { type: 'int?', default: 4, doc: "CPUs to allocate to run rsem"}
  rsem_ram: { type: 'int?', default: 8, doc: "RAM to allocate to run rsem in gb"}
  wf_strand_param: {type: [{type: enum, name: wf_strand_param, symbols: ["default", "rf-stranded", "fr-stranded"]}], doc: "use 'default' for unstranded/auto, 'rf-stranded' if read1 in the fastq read pairs is reverse complement to the transcript, 'fr-stranded' if read1 same sense as transcript"}
  paired: {type: 'boolean?', default: False, doc: "Flag for paired data, separate from wf_strand_param which describes the orientation of paired data [False]"}

outputs:
  matrix_loom: {type: 'File', outputSource: merge_looms/output_file}
  matrix_tar: {type: 'File', outputSource: loom_to_mtx/output_tar}
  alignment_metrics_report: {type: 'File', outputSource: merge_rnaseqc_results/output_file}
  #seurat output

steps:

  strand_parse:
    run: ../tools/expression_parse_strand_param.cwl
    in:
      wf_strand_param: wf_strand_param
    out: [rsem_std, rnaseqc_std, hisat2_std]

  separate_samples:
    run: ../tools/separate_samples.cwl
    in:
      input_dir: input_dir
      paired: paired
    out: [fastq1s, fastq2s, samples]

  build_fastq2_array:
    run: ../tools/build_fastq2_array.cwl
    in:
      fastq1s: separate_samples/fastq1s
      fastq2s: separate_samples/fastq2s
    out: [fastq2_array]

  build_samples_array:
    run: ../tools/build_samples_array.cwl
    in:
      sample_file: separate_samples/samples
    out: [sample_names]

  hisat2_align_genome:
    run: ../tools/hisat2_align.cwl
    scatter: [fastq1, fastq2, output_basename, input_id]
    scatterMethod: dotproduct
    in:
      reference: hisat_genome_ref
      fastq1: separate_samples/fastq1s
      fastq2: build_fastq2_array/fastq2_array
      output_basename: build_samples_array/sample_names
      input_id: build_samples_array/sample_names
      strict:
        valueFrom: ${return Boolean(false)}
      ram: hisat_ram
      cpus: hisat_cpus
    out: [log_file, met_file, bam]

  hisat2_align_trans:
    run: ../tools/hisat2_align.cwl
    scatter: [fastq1, fastq2, output_basename, input_id]
    scatterMethod: dotproduct
    in:
      reference: hisat_trans_ref
      fastq1: separate_samples/fastq1s
      fastq2: build_fastq2_array/fastq2_array
      output_basename: build_samples_array/sample_names
      rna_strandness: strand_parse/hisat2_std
      input_id: build_samples_array/sample_names
      strict:
        valueFrom: ${return Boolean(true)}
      ram: hisat_ram
      cpus: hisat_cpus
    out: [log_file, met_file, bam]

  rnaseqc:
    run: ../tools/rnaseqc.cwl
    scatter: input_bam
    in:
      input_bam: hisat2_align_genome/bam
      collapsed_gtf: rnaseqc_gtf
      strand: strand_parse/rnaseqc_std
      paired: paired
    out: [metrics, gene_TPM, gene_count, exon_count]

  rsem:
    run: ../tools/rsem_calc_express.cwl
    scatter: [input_bam, output_basename]
    scatterMethod: dotproduct
    in:
      input_bam: hisat2_align_trans/bam
      reference: rsem_reference
      output_basename: build_samples_array/sample_names
      paired: paired
      cpus: rsem_cpus
      ram: rsem_ram
    out: [gene_out, isoform_out, cnt_out, model_out, theta_out]

  make_single_loom:
    run: ../tools/make_single_loom.cwl
    scatter: [rsem_matrix, sample_name]
    scatterMethod: dotproduct
    in:
      rsem_matrix: rsem/gene_out
      sample_name: build_samples_array/sample_names
    out: [loom_file]

  merge_rnaseqc_results:
    run: ../tools/merge_rnaseqc_results.cwl
    in:
      output_basename: final_output_basename
      metrics_files: rnaseqc/metrics
    out: [output_file]

  merge_looms:
    run: ../tools/merge_looms.cwl
    in:
      output_basename: final_output_basename
      loom_files: make_single_loom/loom_file
    out: [output_file]

  loom_to_mtx:
    run: ../tools/loom_to_mtx.cwl
    in:
      loom_file: merge_looms/output_file
      output_basename: final_output_basename
    out: [output_tar]

  #Seurat clustering
    #runs on the merged loom file
    #will require even more testing

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4
  - class: 'sbg:AWSInstanceType'
    value: c5.9xlarge
