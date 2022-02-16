<p align="center">
  <img src="docs/kids_first_logo.svg" alt="Kids First repository logo" width="660px" />
</p>
<p align="center">
  <a href="https://github.com/kids-first/kf-template-repo/blob/master/LICENSE"><img src="https://img.shields.io/github/license/kids-first/kf-template-repo.svg?style=for-the-badge"></a>
</p>

# Kids First single cell RNA Pipelines

This repo contains tools and workflows for analyzing single cell RNA (scRNA) data. It includes two workflows one for analyzing 10X data and one for analyzing Smart Seq 2 data.

The repo is currently in alpha phase.

## 10X workflow

The workflow script that runs the tools is `workflows/kf_single_cell_10x.cwl`

The workflow runs [cellranger count](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count),
on fastq files generated by the 10x single cell RNA workflow methodology.
Cell ranger count performs alignment, barcode counting, and filtering.
[SoupX](https://github.com/constantAmateur/SoupX) is used for subtraction of the RNA background
[Scrublet](https://github.com/swolock/scrublet) is used to score and predict doublets
Decontaminated outputs are aggregated using the [Seurat](https://satijalab.org/seurat/) R package from the Satija lab at the New York Genome Center.

### Tools Ran

- Cellranger
- soupX 4.1.0
- scrublet 0.2.3
- Seurat 4.0.4

### Inputs
```yaml
output_basename: {type: string, doc: "basename used to name output files"}
fastq_dirs: {type: 'Directory', doc: "directories of fastqs being run, one from each sample or well"}
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
```

### Outputs
```yaml
count_summary: {type: 'File[]', outputSource: count/output_summary}
bam_out: {type: 'File[]', outputSource: count/bam}
merged_decontam_matrix: {type: 'File', outputSource: merge/merged_matrix}
merged_decontam_object: {type: 'File', outputSource: merge/merged_object}
doublet_histogram: {type: 'File[]', outputSource: scrublet/score_histogram}
```

## Smart Seq 2 Workflow
The workflow uses [HISAT2](http://daehwankimlab.github.io/hisat2/) for alignment, [RNAseQC](https://github.com/getzlab/rnaseqc) to collect sequencing metrics, and [RSEM](https://deweylab.github.io/RSEM/) to calculate gene expression.
The outputs of the workflow are a matrix of cells gene counts in a loom file, a tarball containing a matrix of the gene counts in Matrix Market format, and a tsv with collected sequencing metrics from all cells analyzed.
Basic functionality includes removal of abnormal cells, dimensionality reduction, identification of differentially expressed features, and clustering.

Data are aligned to both the genome and the transcriptome. Genome aligned data are used to collect sequencing metrics and transcriptome aligned data are used to calculate expression.

### Tools Ran
- HISAT2 2.1.0
- RNASeQC 2.4.2
- RSEM 1.3.1

### Inputs
```yaml
final_output_basename: {type: string, doc: "Output basename for workflow output files"}
input_dir: {type: 'Directory', doc: "Directory containing fastq files"}
hisat_genome_ref: {type: 'File', doc: "Hisat 2 genome reference"}
hisat_trans_ref: {type: 'File', doc: "Hisat 2 transcriptome reference"}
rnaseqc_gtf: {type: "File", doc: "gtf file used by RNAseQC", "sbg:suggestedValue": {class: 'File', path: '5d8bb21fe4b0950c4028f852', name: 'gencode.v27.primary_assembly.RNAseQC.gtf'}}
rsem_reference: {type: "File", doc: "RSEM reference file", "sbg:suggestedValue": {class: 'File', path: '5d8bb21fe4b0950c4028f851', name: 'RSEM_GENCODE27.tar.gz'}}
cpus: { type: 'int?', default: 4, doc: "CPUs to allocate to call task"}
ram: { type: 'int?', default: 8, doc: "RAM to allocate to call task in gb"}
wf_strand_param: {type: [{type: enum, name: wf_strand_param, symbols: ["default", "rf-stranded", "fr-stranded"]}], doc: "use 'default' for unstranded/auto, 'rf-stranded' if read1 in the fastq read pairs is reverse complement to the transcript, 'fr-stranded' if read1 same sense as transcript"}
paired: {type: 'boolean?', default: False, doc: "Flag for paired data, separate from wf_strand_param which describes the orientation of paired data [False]"}
```

### Outputs
```yaml
matrix_loom: {type: 'File', outputSource: merge_looms/output_file}
matrix_tar: {type: 'File', outputSource: loom_to_mtx/output_tar}
alignment_metrics_report: {type: 'File', outputSource: merge_rnaseqc_results/output_file}
```
