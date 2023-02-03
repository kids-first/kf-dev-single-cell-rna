## Smart Seq 2 Workflow
The workflow uses [HISAT2](http://daehwankimlab.github.io/hisat2/) for alignment, [RNAseQC](https://github.com/getzlab/rnaseqc) to collect sequencing metrics, and [RSEM](https://deweylab.github.io/RSEM/) to calculate gene expression.
The HISAT2 + RSEM run paramaters were guided by the [Human Cell Atlas Data Portal](https://data.humancellatlas.org/pipelines/smart-seq2-workflow).
The outputs of the workflow are a matrix of cells gene counts in a loom file, a tarball containing a matrix of the gene counts in Matrix Market format, and a tsv with collected sequencing metrics from all cells analyzed.
Basic functionality includes removal of abnormal cells, dimensionality reduction, identification of differentially expressed features, and clustering.

Data are aligned to both the genome and the transcriptome. Genome aligned data are used to collect sequencing metrics and transcriptome aligned data are used to calculate expression.

### Tools Ran
- HISAT2 2.1.0
- RNASeQC 2.4.2
- RSEM 1.3.1

### Inputs
 - final_output_basename: Output basename for workflow output files
 - input_dir: Directory containing fastq files
 - hisat_genome_ref: Hisat 2 genome reference
 - hisat_trans_ref: Hisat 2 transcriptome reference
 - rnaseqc_gtf: gtf file used by RNAseQC, recommend gencode.v39.primary_assembly.rnaseqc.stranded.gtf
 - rsem_reference: RSEM reference file, recommend RSEM_GENCODE39.tar.gz
 - cpus: CPUs to allocate to call task
 - ram: RAM to allocate to call task in gb
 - wf_strand_param: use 'default' for unstranded/auto, 'rf-stranded' if read1 in the fastq read pairs is reverse complement to the transcript, 'fr-stranded' if read1 same sense as transcript
 - paired: Flag for paired data, separate from wf_strand_param which describes the orientation of paired data [False]

### Outputs
```yaml
matrix_loom: {type: 'File', outputSource: merge_looms/output_file}
matrix_tar: {type: 'File', outputSource: loom_to_mtx/output_tar}
alignment_metrics_report: {type: 'File', outputSource: merge_rnaseqc_results/output_file}
```

### Custom reference building
#### hisat_genome_ref
The name is confusing as it includes a transcript refrence as part of th ebuild process, but uses the fasta genome fa. Built using tools/hisat2_build_index.cwl, with the following preprocessing steps to generate inputs
 - Used GENCODE39 primary reference and tools/hisat2_format_gtf_ref.cwl tool to generate `exon` and `splice` inputs
 - Followed instructions from https://daehwankimlab.github.io/hisat2/howto/#build-hgfm-index-with-snps, but skipped awk command as 
    used snp151 instead of 144 to generate `snp` input locally
#### hisat_trans_ref
Used for alignment against a transcriptome fasta only for RSEM quant.
Built using tools/hisat2_build_index.cwl with the fasta file inside the RSEM tar ball reference used as fasta input. No other inputs aside from output basename given.