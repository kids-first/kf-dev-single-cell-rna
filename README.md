<p align="center">
  <img src="docs/kids_first_logo.svg" alt="Kids First repository logo" width="660px" />
</p>
<p align="center">
  <a href="https://github.com/kids-first/kf-template-repo/blob/master/LICENSE"><img src="https://img.shields.io/github/license/kids-first/kf-template-repo.svg?style=for-the-badge"></a>
</p>

# Kids First single cell RNA Pipelines

This repo contains tools and workflows for analyzing single cell RNA (scRNA) data from 10X and SmartSeq2

## Production Workflows
These workflows are considered acceptable for Single Cell Anaylsis
### [10X STAR Solo 2.7.10b Alignment Workflow](docs/10X_STAR_Solo_alignment.md)
Aligns and quantifies 10X single cell data using STAR Solo

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
