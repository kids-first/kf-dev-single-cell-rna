# KF Cellranger Count/Multi Workflow
This workflow takes 10X FASTQ files and processes them either through `count` or `multi` depending on the `mode` flag.

## Inputs
### General Required
 - `mode`: `count` or `multi`. In general, `count` is good for single sample, `multi` for more than one sample, required for Flex chemistries.
 See [Choosing a pipeline](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-choosing-a-pipeline) for more details.
 - `reads`: R1 FASTQ files
 - `mates`: R2 FASTQ files
 - `transcriptome_dir`: Unarchived cell ranger reference if `transriptome_tar` not given
 - `transcriptome_tar`: Archived cell ranger reference if `transriptome_dir` not given. See [references](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads) for more info
 - `create_bam`: Flag to create alignment bams along side counts files.
 Default `false` as it's slower, but may be useful as a demultiplexing mechanism.
### Count-specific required
 - `sample`: Sample name to use
### Count-specific if files exist
 - `indices`: `I1` FASTQ files
### Multi-specific required
 - `sample_sheet`: CSV file with columns: sample_id,probe_barcode_ids
This is to fill out the [sample] section of the multi config.

 - `library_fastq_id`: Should be the prefix of the input FASTQ up until the `_S#_L##_R#_001.fastq.gz` part
 - `feature_types`: Default is `Gene Expression`. See [docs](https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-multi-config-csv-opts#feature-types)
### Multi-specific Flex
See [Multi Config](https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-multi-config-csv-opts) for details
- `probe_set`: Probe set file obtainable from `references`
