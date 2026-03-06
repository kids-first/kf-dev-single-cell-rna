<p align="center">
  <img src="docs/kids_first_logo.svg" alt="Kids First repository logo" width="660px" />
</p>
<p align="center">
  <a href="https://github.com/kids-first/kf-dev-single-cell-rna/blob/master/LICENSE"><img src="https://img.shields.io/github/license/kids-first/kf-dev-single-cell-rna.svg?style=for-the-badge"></a>
</p>

# Kids First single cell RNA Pipelines

This repo contains tools and workflows for analyzing single cell RNA (scRNA) data from 10X and SmartSeq2

## Production Workflows
These workflows are considered acceptable for Single Cell Analysis
### [10X STAR Solo 2.7.10b Alignment Workflow](docs/10X_STAR_Solo_alignment.md)
Aligns and quantifies 10X single cell data using STAR Solo

### [KF Cellranger Count/Multi Workflow](docs/CELLRANGER_COUNT_MULTI.md)
Nextflow workflow to generate gene expression (gex) counts files using `count` or `multi` modes.
Can also generate alignments as well as cell ranger cluster output.

## Beta Phase Workflows
These workflows are currently in beta phase

### [SWANS 2.0: Single-entity Workflow ANalysiS Pipeline](docs/SWANS_PRELIMINARY_ANALYSIS.md)
Complete QC and analysis pipeline. https://doi.org/10.1101/2025.05.14.654073. It is a three part process consisting of the following:
- [Data Cleanup and QC](docs/SWANS_DATA_CLEANUP_AND_QC.md)
  - Run tools like soupX and doubletFinder, Seurat filtering, etc. May need to re-run to fine-tune after review
  - Both [nextflow](nf/SWANS_DATA_QC/main.nf) and [CAVATICA-ready](nf/SWANS_DATA_QC/sb_nextflow_schema.yaml)
- [Interactive Report](docs/SWANS_INTERACTIVE_REPORT.md)
  - Once QC params tuned, generate preliminary results for study and review
  - Both [nextflow](nf/SWANS_INTERACTIVE/main.nf) and [CAVATICA-ready](nf/SWANS_INTERACTIVE/sb_nextflow_schema.yaml)
- [Post Analysis](docs/SWANS_POST_ANALYSIS.md)
  - aka, "Final Analysis". After careful review of interactive report results, zero in on specific params
  - [nextflow](nf/SWANS_INTERACTIVE/main.nf) in alpha

## Alpha Phase Workflows
These workflows are currently in alpha phase

### [10X Cell Ranger 6.2.1 Alignment Workflow](docs/10X_cell_ranger_alignment.md)
Aligns and quantifies 10X single cell data using Cell Ranger
### [10X Refinement Workflow](docs/10X_refinement.md)
Filters and cleans up results from both Cell Ranger and STAR Solo outputs

### [Smart Seq 2 Workflow](docs/SMART_SEQ2.md)
The workflow uses [HISAT2](http://daehwankimlab.github.io/hisat2/) for alignment, [RNAseQC](https://github.com/getzlab/rnaseqc) to collect sequencing metrics, and [RSEM](https://deweylab.github.io/RSEM/) to calculate gene expression.
The outputs of the workflow are a matrix of cells gene counts in a loom file, a tarball containing a matrix of the gene counts in Matrix Market format, and a tsv with collected sequencing metrics from all cells analyzed.
Basic functionality includes removal of abnormal cells, dimensionality reduction, identification of differentially expressed features, and clustering.

Data are aligned to both the genome and the transcriptome. Genome aligned data are used to collect sequencing metrics and transcriptome aligned data are used to calculate expression.
