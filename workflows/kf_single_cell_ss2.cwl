cwlVersion: v1.0
class: Workflow
id: kf_singe_cell_ss2
label: Kids First DRC single cell RNA Smart Seq 2 Workflow
doc: |-
  # KFDRC single cell RNA Smart Seq 2 workflow
  Workflow for aligning and counting single cell RNA reads generated by Smart Seq 2.

  ![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

inputs:
  final_output_basename: {type: string, doc: "Output basename for workflow output files"}
  fastq1s: {type: 'File[]', doc: "Array of fastq 1s to align"}
  fastq2s: {type: ['null', 'File[]?'], doc: "Array of fastq 2s to align"}
  sample_names: {type: 'string[]', doc: "Array of sample names"}
  hisat_genome_ref: {type: 'File', doc: "Hisat 2 genome reference"}
  hisat_trans_ref: {type: 'File', doc: "Hisat 2 transcriptome reference"}
  rnaseqc_gtf: {type: 'File', doc: "gtf file used by RNAseQC", sbg:suggestedValue: {class: 'File', path: '5d8bb21fe4b0950c4028f852', name: 'gencode.v27.primary_assembly.RNAseQC.gtf'}}
  rsem_reference: {type: 'File', doc: "RSEM reference file", sbg:suggestedValue: {class: 'File', path: '5d8bb21fe4b0950c4028f851', name: 'RSEM_GENCODE27.tar.gz'}}
  cpus: { type: 'int?', default: 4, doc: "CPUs to allocate to call task"}
  ram: { type: 'int?', default: 8, doc: "RAM to allocate to call task in gb"}
  wf_strand_param: {type: [{type: enum, name: wf_strand_param, symbols: ["default", "rf-stranded", "fr-stranded"]}], doc: "use 'default' for unstranded/auto, 'rf-stranded' if read1 in the fastq read pairs is reverse complement to the transcript, 'fr-stranded' if read1 same sense as transcript"}

outputs:
  matrix_loom: {type: 'File', outputSource: merge_looms/output_file}
  alignment_metrics_report: {type: 'File', outputSource: merge_rnaseqc_results/output_file}
  #seurat output

steps:

  build_fastq2_array:
    run: ../tools/build_fastq2_array.cwl
    in:
      fastq1s: fastq1s
      fastq2s: fastq2s
    out: [fastq2s]

  strand_parse:
    run: ../tools/expression_parse_strand_param.cwl
    in:
      wf_strand_param: wf_strand_param
    out: [rsem_std, rnaseqc_std]

  hisat2_align_genome:
    run: ../tools/hisat2_align.cwl
    scatter: [fastq1, fastq2, output_basename, input_id]
    scatterMethod: dotproduct
    in:
      reference: hisat_genome_ref
      fastq1: fastq1s
      fastq2: build_fastq2_array/fastq2s
      output_basename: sample_names
      input_id: sample_names
      strict:
        valueFrom: ${return Boolean(false)}
      ram: ram
      cpus: cpus
    out: [log_file, met_file, bam]

  hisat2_align_trans:
    run: ../tools/hisat2_align.cwl
    scatter: [fastq1, fastq2, output_basename, input_id]
    scatterMethod: dotproduct
    in:
      reference: hisat_trans_ref
      fastq1: fastq1s
      fastq2: build_fastq2_array/fastq2s
      output_basename: sample_names
      input_id: sample_names
      strict:
        valueFrom: ${return Boolean(true)}
      ram: ram
      cpus: cpus
    out: [log_file, met_file, bam]

  rnaseqc:
    run: ../tools/rnaseqc.cwl
    scatter: input_bam
    in:
      input_bam: hisat2_align_genome/bam
      collapsed_gtf: rnaseqc_gtf
      strand: strand_parse/rnaseqc_std
      paired:
        valueFrom: ${
          if (inputs.fastq2s) { return Boolean(true)}
          else { return Boolean(false)}
          }
      ram: ram
      cpus: cpus
    out: [metrics, gene_TPM, gene_count, exon_count]

  rsem:
    run: ../tools/rsem_calc_express.cwl
    scatter: [input_bam, output_basename]
    scatterMethod: dotproduct
    in:
      input_bam: hisat2_align_trans/bam
      reference: rsem_reference
      output_basename: sample_names
      paired:
        valueFrom: ${
          if (inputs.fastq2s) { return Boolean(true)}
          else { return Boolean(false)}
          }
    out: [gene_out, isoform_out, cnt_out, model_out, theta_out]

  make_single_loom:
    run: ../tools/make_single_loom.cwl
    scatter: [rsem_matrix, sample_name]
    scatterMethod: dotproduct
    in:
      rsem_matrix: rsem/gene_out
      sample_name: sample_names
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

  #Seurat clustering
    #runs on the merged loom file
    #will require even more testing

$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 5
