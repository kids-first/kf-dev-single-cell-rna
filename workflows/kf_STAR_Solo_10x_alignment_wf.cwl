cwlVersion: v1.2
class: Workflow
id: kf-star-solo-10x-align-wf
label: "KFDRC Single Cell RNA STAR Solo 10x Alignment Workflow"
doc: |
  # 10X Alignment Workflow

  The workflow script that runs the tools is `workflows/kf_STAR_Solo_10x_alignment.cwl`

  The workflow runs [STAR Solo](https://github.com/alexdobin/STAR/blob/2.7.10b/docs/STARsolo.md),
  on fastq files generated by the 10x single cell RNA workflow methodology.
  Star Solo performs alignment, barcode counting, and filtering.
  Quote from the STAR repo manual pdf:
  >STARsolo output is designed to be a drop-in replacement for 10X CellRanger gene quantification output. It follows CellRanger logic for cell barcode whitelisting and UMI deduplication, and produces nearly identical gene counts in the same format. At the same time STARsolo is 10 times faster than the CellRanger.

  Some custom tools are used to create Cell Ranger-like output file structure and h5 files for downstream refinement compatibility
  Output QC is based on [this tutorial](https://github.com/hbctraining/scRNA-seq_online/blob/scRNAseq/lessons/04_SC_quality_control.md) and [this doc](./10X_SEURAT_HBC_SCRNA_QC.md) summarizes the process and outputs
  ## Software

  - STAR Solo 2.7.10b
  - Seurat 4.3.0.1

  ## Inputs
  ### multi-step
   - `output_basename`: basename used to name output files
   - `sample_name`: used as prefix for finding fastqs to analyze, e.g. 1k_PBMCs_TotalSeq_B_3p_LT_antibody if the names of the underlying fastqs are of the form 1k_PBMCs_TotalSeq_B_3p_LT_antibody_S1_L001_I1_001.fastq.gz, one per input fastq in the same order
  ### STAR Solo
   - `outSAMattrRGline`: Set if outputting bam, with TABS SEPARATING THE TAGS, format is: ID:sample_name LB:aliquot_id PL:platform SM:BSID for example ID:7316-242	LB:750189	PL:ILLUMINA	SM:BS_W72364MN
   - `genomeDir`: Tar gzipped reference that will be unzipped at run time
   - `readFilesIn1`: Input fastq file(s), gzipped or uncompressed
   - `readFilesIn2` R2 or 'mates' reads file(s), gzipped or uncompressed
   - `runThreadN` Num threads for STAR Solo to use, default: 16
   - `solo_type`: type of single-cell RNAseq
     - CB_UMI_Simple: (a.k.a. Droplet) one UMI and one Cell Barcode of fixed length in read2, e.g. Drop-seq and 10X Chromium
     - CB UMI_Complex: multiple Cell Barcodes of varying length, one UMI of fixed length and one adapter sequence of fixed length are allowed in read2 only (e.g. inDrop, ddSeq)
     - CB samTagOut: output Cell Barcode as CR and/or CB SAm tag. No UMI counting. –readFilesIn cDNA read1 [cDNA read2 if paired-end] CellBarcode read . Requires -outSAMtype BAM Unsorted [and/or SortedByCoordinate]
     - Smart-seq: each cell in a separate FASTQ (paired- or single-end), barcodes are corresponding read-groups, no UMI sequences, alignments deduplicated according to alignment start and end (after extending soft-clipped bases)
     - default: "CB_UMI_Simple"
   - `soloCBwhitelist`: file with whitelist of cell barcodes
   - `soloUMIlen`: UMI length, default: 12
   - `clipAdapterType`: adapter clipping type.
     - Hamming: adapter clipping based on Hamming distance, with the number of mismatches controlled by -clip5pAdapterMMp
     - CellRanger4: 5p and 3p adapter clipping similar to CellRanger4. Utilizes Opalpackage by Martin Sosic: https://github.com/Martinsos/opal
     - None: no adapter clipping, all other clip* parameters are disregarded
     - default: "CellRanger4"
   - `outFilterScoreMin`: alignment will be output only if its score is higher than or equal to this value, default: 30
   - `soloCBmatchWLtype`: matching the Cell Barcodes to the WhiteList.
     - Exact: only exact matches allowed
     - 1MM: only one match in whitelist with 1 mismatched base allowed. Allowed CBs have to have at least one read with exact match.
     - 1MM_multi: multiple matches in whitelist with 1 mismatched base allowed, posterior probability calculation is used choose one of the matches. Allowed CBs have to have at least one read with exact match. This option matches best with CellRanger 2.2.0
     - 1MM_multi_pseudocounts: same as 1MM Multi, but pseudocounts of 1 are added to all whitelist barcodes
     - 1MM_multi_Nbase_pseudocounts: same as 1MM multi pseudocounts, multimatching to WL is allowed for CBs with N-bases. This option matches best with CellRanger >= 3.0.0
     - EditDist_2: allow up to edit distance of 3 fpr each of the barcodes. May include one deletion + one insertion. Only works with -soloType CB UMI Complex. Matches to multiple passlist barcdoes are not allowed. Similar to ParseBio Split-seq pipeline
     - default: "1MM_multi_Nbase_pseudocounts"
   - `soloUMIfiltering`: type of UMI filtering (for reads uniquely mapping to genes)
     - basic filtering: remove UMIs with N and homopolymers (similar to CellRanger 2.2.0)
     - MultiGeneUMI: basic + remove lower-count UMIs that map to more than one gene
     - MultiGeneUMI_All: basic + remove all UMIs that map to more than one gene
     - MultiGeneUMI_CR: basic + remove lower-count UMIs that map to more than one gene, matching CellRanger > 3.0.0. Only works with -soloUMIdedup 1MM_CR
     - default: "MultiGeneUMI_CR"
   - `soloUMIdedup`: type of UMI deduplication (collapsing) algorithm)
     - 1MM_All: all UMIs with 1 mismatch distance to each other are collapsed (i.e. counted once)
     - 1MM_Directional_UMItools: follows the ”directional” method from the UMI-tools by Smith, Heger and Sudbery (Genome Research 2017)
     - 1MM_Directional: same as 1MM Directional UMItools, but with more stringent criteria for duplicate UMIs
     - Exact: only exactly matching UMIs are collapsed
     - NoDedup: no deduplication of UMIs, count all reads
     - 1MM_CR: CellRanger2-4 algorithm for 1MM UMI collapsing
     - default: "1MM_CR"
   - `soloFeatures`: genomic features for which the UMI counts per Cell Barcode are collected. **Recommend for 10X at least `Gene` for whole cell input, `GeneFull` for snRNA input**
     - Gene: genes: reads match the gene transcript
     - SJ: splice junctions: reported in SJ.out.tab
     - GeneFull: full gene (pre-mRNA): full gene (pre-mRNA): count all reads overlapping genes' exons and introns
     - GeneFull_ExonOverIntron: count all reads overlapping genes' exons and introns: prioritize 100% overlap with exons
     - GeneFull_Ex50pAS: full gene (pre-RNA): count all reads overlapping genes' exons and introns: prioritize >50% overlap with exons. Do not count reads with 100% exonic overlap in the antisense direction
     - default: "GeneFull" 
   - `soloCellFilter`: cell filtering type and parameters
     - None: do not output filtered cells
     - TopCells: only report top cells by UMI count, followed by the exact number of cells
     - CellRanger2.2: simple filtering of CellRanger 2.2. Can be followed by numbers: number of expected cells, robust maximum percentile for UMI count, maximum to minimum ratio for UMI count. The harcoded values are from CellRanger: nExpectedCells=3000; maxPercentile=0.99; maxMinRatio=10
     - EmptyDrops_CR: EmptyDrops filtering in CellRanger flavor. Please cite the original EmptyDrops paper: A.T.L Lun et al, Genome Biology, 20, 63 (2019): https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y, Can be followed by 10 numeric parameters: nExpectedCells maxPercentile maxMinRatio indMin indMax umiMin umiMinFracMedian candMaxN FDR simN. The harcoded values are from CellRanger: 3000 0.99 10 45000 90000 500 0.01 20000 0.01 10000
     - default: "EmptyDrops_CR"
    outSAMtype: type of SAM/BAM output. None: no SAM/BAM output. Otherwise, first word is output type (BAM or SAM), second is sort type (Unsorted or SortedByCoordinate)
     - default: "None"
  ### seurat hbc qc
   - `qc_min_umi`: minimum number of umi for cell-level filtering
   - `qc_min_genes`: minimum number of genes for cell-level filtering
   - `qc_min_complexity`: minimum novelty score (log10GenesPerUMI)
   - `qc_max_mito_ratio`: maximum ratio mitochondrial reads per cell
   - `qc_min_gene_prevalence`: Minimum number of cells a gene must be expressed in to keep after filtering

  ## Outputs
   - `star_solo_counts_dir`: Tar gzipped counts output from STAR Solo
   - `star_solo_bam`: If flag given, aligned reads file
   - `star_solo_matrix_filtered`: Filtered counts matrix in Cell Ranger-style h5 format
   - `star_solo_matrix_raw`: Raw counts matrix in Cell Ranger-style h5 format
   - `star_solo_log_final_out`: STAR align summary stats
   - `star_solo_junctions`: STAR splice junction result file
   - `star_solo_cr_mimic_counts`: Tar ball of Cell Ranger-style counts dir from STAR Solo
   - `qc_plots`: Pre and post filtering metrics PDF plots
   - `qc_boxplot_stats`: Pre and post filtering boxplot stats TSV
   - `cell_counts`: Pre and post filtering cell counts TSV
   - `seurat_prefilter_data`: Seurat Rdata object with prefilter counts and metrics
   - `seurat_filtered_data`: Seurat Rdata object with basic filter counts and metrics
   - `variable_features_plot`: PDF with a dot plot of variable genes with top 15 labeled
  ## Appendix: Seurat HBC QC Output
  QC Outputs are based primarily on the training materials provided by the [Harvard Chan Bioinformatics Core](https://github.com/hbctraining/scRNA-seq_online/blob/scRNAseq/lessons/04_SC_quality_control.md).
  An overview of how the QC was performed and overview of the outputs are provided [here](./10X_SEURAT_HBC_SCRNA_QC.md)
requirements:
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
inputs:
  # multi-step
  output_basename: {type: string, doc: "basename used to name output files"}
  sample_name: {type: string, doc: "used as prefix for labeling data for downstream anaylsis"}
  genomeDir: {type: File, doc: "Tar gzipped reference that will be unzipped at run time"}
  readFilesIn1: {type: 'File[]', doc: "Input fastq file(s), gzipped or uncompressed"}
  readFilesIn2: {type: 'File[]', doc: "R2 or 'mates' reads file(s), gzipped or uncompressed"}
  runThreadN: {type: 'int?', doc: "Num threads for STAR Solo to use", default: 16}
  outSAMattrRGline: {type: 'string?', doc: "Set if outputting bam, with TABS SEPARATING THE TAGS, format is: ID:sample_name LB:aliquot_id
      PL:platform SM:BSID for example ID:7316-242 LB:750189 PL:ILLUMINA SM:BS_W72364MN"}
  twopassMode: {type: ['null', {type: enum, name: twopassMode, symbols: ["Basic", "None"]}], default: "Basic", doc: "Enable two pass
      mode to detect novel splice events. Default is Basic (on)."}
  solo_type: {type: ['null', {type: enum, name: soloType, symbols: ["CB_UMI_Simple", "CB UMI Complex", "CB samTagOut", "SmartSeq"]}],
    doc: "type of single-cell RNAseq. CB_UMI_Simple: (a.k.a. Droplet) one UMI and one Cell Barcode of fixed length in read2, e.g.
      Drop-seq and 10X Chromium, CB UMI_Complex: multiple Cell Barcodes of varying length, one UMI of fixed length and one adapter
      sequence of fixed length are allowed in read2 only (e.g. inDrop, ddSeq), CB samTagOut: output Cell Barcode as CR and/or CB SAm
      tag. No UMI counting. –readFilesIn cDNA read1 [cDNA read2 if paired-end] CellBarcode read . Requires -outSAMtype BAM Unsorted
      [and/or SortedByCoordinate] Smart-seq: each cell in a separate FASTQ (paired- or single-end), barcodes are corresponding read-groups,
      no UMI sequences, alignments deduplicated according to alignment start and end (after extending soft-clipped bases)", default: "CB_UMI_Simple"}
  soloCBwhitelist: {type: File, doc: "file with whitelist of cell barcodes"}
  soloUMIlen: {type: 'int?', doc: "UMI length", default: 12}
  clipAdapterType: {type: ['null', {type: enum, name: clipAdapterType, symbols: ["Hamming", "CellRanger4", "None"]}], doc: "adapter
      clipping type. Hamming: adapter clipping based on Hamming distance, with the number of mismatches controlled by -clip5pAdapterMMp.
      CellRanger4: 5p and 3p adapter clipping similar to CellRanger4. Utilizes Opalpackage by Martin Sosic: https://github.com/Martinsos/opal,None:
      no adapter clipping, all other clip* parameters are disregarded", default: "CellRanger4"}
  outFilterScoreMin: {type: 'int?', doc: "alignment will be output only if its score is higher than or equal to this value", default: 30}
  soloCBmatchWLtype: {type: ['null', {type: enum, name: clipAdapterType, symbols: ["Exact", "1MM", "1MM_multi", "1MM_multi_pseudocounts",
          "1MM_multi_Nbase_pseudocounts", "EditDist_2"]}], doc: "matching the Cell Barcodes to the WhiteList. Exact: only exact matches
      allowed. 1MM: only one match in whitelist with 1 mismatched base allowed. Allowed CBs have to have at least one read with exact
      match. 1MM_multi: multiple matches in whitelist with 1 mismatched base allowed, posterior probability calculation is used choose
      one of the matches. Allowed CBs have to have at least one read with exact match. This option matches best with CellRanger 2.2.0.
      1MM_multi_pseudocounts: same as 1MM Multi, but pseudocounts of 1 are added to all whitelist barcodes. 1MM_multi_Nbase_pseudocounts:
      same as 1MM multi pseudocounts, multimatching to WL is allowed for CBs with N-bases. This option matches best with CellRanger
      >= 3.0.0 EditDist_2: allow up to edit distance of 3 fpr each of the barcodes. May include one deletion + one insertion. Only
      works with -soloType CB UMI Complex. Matches to multiple passlist barcdoes are not allowed. Similar to ParseBio Split-seq pipeline.",
    default: "1MM_multi_Nbase_pseudocounts"}
  soloUMIfiltering: {type: ['null', {type: enum, name: soloUMIfiltering, symbols: ["-", "MultiGeneUMI", "MultiGeneUMI_All", "MultiGeneUMI_CR"]}],
    doc: "type of UMI filtering (for reads uniquely mapping to genes). -: basic filtering: remove UMIs with N and homopolymers (similar
      to CellRanger 2.2.0). MultiGeneUMI: basic + remove lower-count UMIs that map to more than one gene. MultiGeneUMI_All: basic
      + remove all UMIs that map to more than one gene. MultiGeneUMI_CR: basic + remove lower-count UMIs that map to more than one
      gene, matching CellRanger > 3.0.0. Only works with -soloUMIdedup 1MM_CR", default: "MultiGeneUMI_CR"}
  soloUMIdedup: {type: ['null', {type: enum, name: soloUMIdedup, symbols: ["1MM_All", "1MM_Directional_UMItools", "1MM_Directional",
          "Exact", "NoDedup", "1MM_CR"]}], doc: "type of UMI deduplication (collapsing) algorithm). 1MM_All: all UMIs with 1 mismatch
      distance to each other are collapsed (i.e. counted once). 1MM_Directional_UMItools: follows the 'directional' method from the
      UMI-tools by Smith, Heger and Sudbery (Genome Research 2017). 1MM_Directional: same as 1MM Directional UMItools, but with more
      stringent criteria for duplicate UMIs. Exact: only exactly matching UMIs are collapsed NoDedup: no deduplication of UMIs, count
      all reads. 1MM_CR: CellRanger2-4 algorithm for 1MM UMI collapsing.", default: "1MM_CR"}
  soloFeatures: {type: ['null', {type: enum, name: soloFeatures, symbols: ["Gene", "SJ", "GeneFull", "GeneFull_ExonOverIntron", "GeneFull_Ex50pAS"]}],
    doc: "genomic features for which the UMI counts per Cell Barcode are collected. Gene: genes: reads match the gene transcript.
      SJ: splice junctions: reported in SJ.out.tab. GeneFull: full gene (pre-mRNA): full gene (pre-mRNA): count all reads overlapping
      genes' exons and introns. GeneFull_ExonOverIntron: count all reads overlapping genes' exons and introns: prioritize 100% overlap
      with exons. GeneFull_Ex50pAS: full gene (pre-RNA): count all reads overlapping genes' exons and introns: prioritize >50% overlap
      with exons. Do not count reads with 100% exonic overlap in the antisense direction.", default: "GeneFull"}
  soloCellFilter: {type: ['null', {type: enum, name: soloCellFilter, symbols: ["None", "TopCells", "CellRanger2.2", "EmptyDrops_CR"]}],
    doc: "cell filtering type and parameters. None: do not output filtered cells. TopCells: only report top cells by UMI count, followed
      by the exact number of cells. CellRanger2.2: simple filtering of CellRanger 2.2. Can be followed by numbers: number of expected
      cells, robust maximum percentile for UMI count, maximum to minimum ratio for UMI count. The harcoded values are from CellRanger:
      nExpectedCells=3000; maxPercentile=0.99; maxMinRatio=10. EmptyDrops_CR: EmptyDrops filtering in CellRanger flavor. Please cite
      the original EmptyDrops paper: A.T.L Lun et al, Genome Biology, 20, 63 (2019): https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y,
      Can be followed by 10 numeric parameters: nExpectedCells maxPercentile maxMinRatio indMin indMax umiMin umiMinFracMedian candMaxN
      FDR simN. The harcoded values are from CellRanger: 3000 0.99 10 45000 90000 500 0.01 20000 0.01 10000", default: "EmptyDrops_CR"}
  soloMultiMappers: {type: 'string[]?', doc: "Possible one or more values: 'Unique', 'Uniform', 'PropUnique', 'EM', 'Rescue'. Including
      multi-gene reads allows for more accurate gene quantification and, more importantly, enables detection of gene expression from
      certain classes of genes that are supported only by multi-gene reads, such as overlapping genes and highly similar paralog families.
      Unique: software default, count only reads that map to unique genes Uniform: uniformly distributes the multi-gene UMIs to all
      genes in its gene set. Each gene gets a fractional count of 1/N_genes, where N_genes is the number of genes in the set. This
      is the simplest possible option, and it offers higher sensitivity for gene detection at the expense of lower precision PropUnique:
      distributes the multi-gene UMIs proportionally to the number of unique UMIs per gene. UMIs that map to genes that are not supported
      by unique UMIs are distributed uniformly EM: uses Maximum Likelihood Estimation (MLE) to distribute multi-gene UMIs among their
      genes, taking into account other UMIs (both unique- and multi-gene) from the same cell (i.e. with the same CB). Expectation-Maximization
      (EM) algorithm is used to find the gene expression values that maximize the likelihood function. Recovering multi-gene reads
      via MLE-EM model was previously used to quantify transposable elements in bulk RNA-seq {TEtranscripts} and in scRNA-seq {Alevin;
      Kallisto-bustools}. Rescue: distributes multi-gene UMIs to their gene set proportionally to the sum of the number of unique-gene
      UMIs and uniformly distributed multi-gene UMIs in each gene Mortazavi et al. It can be thought of as the first step of the EM
      algorithm"}
  raw_count_choice: {type: ['null', {type: enum, name: raw_count_choice, symbols: ["Unique", "Uniform", "PropUnique", "EM", "Rescue"]}],
    doc: "Based on `soloMultiMappers`, if you wish to include/handle multi-gene hits in downstream anaylsis instead of default (ignore
      multi-gene mappers), pick the method you want to use", default: "Unique"}
  outSAMtype: {type: ['null', {type: enum, name: outSAMtype, symbols: ["BAM Unsorted", "None", "BAM SortedByCoordinate", "SAM Unsorted",
          "SAM SortedByCoordinate"]}], default: "None", doc: "type of SAM/BAM output. None: no SAM/BAM output. Otherwise, first word
      is output type (BAM or SAM), second is sort type (Unsorted or SortedByCoordinate)"}
  # Seurat HBC QC
  qc_min_umi: {type: 'int?', doc: "minimum number of umi for cell-level filtering", default: 500}
  qc_min_genes: {type: 'int?', doc: "minimum number of genes for cell-level filtering", default: 250}
  qc_min_complexity: {type: 'float?', doc: "minimum novelty score (log10GenesPerUMI)", default: 0.8}
  qc_max_mito_ratio: {type: 'float?', doc: "maximum ratio mitochondrial reads per cell", default: 0.2}
  qc_min_gene_prevalence: {type: 'int?', doc: "Minimum number of cells a gene must be expressed in to keep after filtering", default: 10}
outputs:
  star_solo_counts_dir: {type: File, outputSource: tar_solo_count_outdir/output, doc: "Tar gzipped counts output from STAR Solo"}
  star_solo_bam: {type: 'File?', outputSource: star_solo_align/genomic_bam_out, doc: "If flag given, aligned reads file"}
  star_solo_matrix_filtered: {type: File, outputSource: create_h5_output/filtered_converted_h5, doc: "Filtered counts matrix in Cell
      Ranger-style h5 format"}
  star_solo_matrix_raw: {type: File, outputSource: create_h5_output/raw_converted_h5, doc: "Raw counts matrix in Cell Ranger-style
      h5 format"}
  star_solo_log_final_out: {type: File, outputSource: star_solo_align/log_final_out, doc: "STAR align summary stats"}
  star_solo_junctions: {type: File, outputSource: star_solo_align/junctions_out, doc: "STAR splice junction result file"}
  star_solo_cr_mimic_counts: {type: File, outputSource: tar_solo_cr_mimic_dir/output, doc: "Tar ball of Cell Ranger-style counts dir
      from STAR Solo"}
  qc_plots: {type: File, outputSource: seurat_hbc_qc/qc_plots, doc: "Pre and post filtering metrics PDF plots"}
  qc_boxplot_stats: {type: File, outputSource: seurat_hbc_qc/qc_boxplot_stats, doc: "Pre and post filtering boxplot stats TSV"}
  cell_counts: {type: File, outputSource: seurat_hbc_qc/cell_counts, doc: "Pre and post filtering cell counts TSV"}
  seurat_prefilter_data: {type: File, outputSource: seurat_hbc_qc/seurat_prefilter_data, doc: "Seurat Rdata object with prefilter
      counts and metrics"}
  seurat_filtered_data: {type: File, outputSource: seurat_hbc_qc/seurat_filtered_data, doc: "Seurat Rdata object with basic filter
      counts and metrics"}
  variable_features_plot: {type: File, outputSource: seurat_hbc_qc/variable_features_plot, doc: "PDF with a dot plot of variable genes
      with top 15 labeled"}
steps:
  star_solo_align:
    run: ../tools/star_solo_2.7.10b.cwl
    in:
      outSAMattrRGline: outSAMattrRGline
      twopassMode: twopassMode
      genomeDir: genomeDir
      readFilesIn1: readFilesIn1
      readFilesIn2: readFilesIn2
      runThreadN: runThreadN
      outFileNamePrefix: output_basename
      solo_type: solo_type
      soloCBwhitelist: soloCBwhitelist
      soloUMIlen: soloUMIlen
      clipAdapterType: clipAdapterType
      outFilterScoreMin: outFilterScoreMin
      soloCBmatchWLtype: soloCBmatchWLtype
      soloUMIfiltering: soloUMIfiltering
      soloUMIdedup: soloUMIdedup
      soloFeatures: soloFeatures
      soloCellFilter: soloCellFilter
      soloMultiMappers: soloMultiMappers
      outSAMtype: outSAMtype
    out: [log_progress_out, log_out, log_final_out, genomic_bam_out, junctions_out, counts_dir]
  create_h5_output:
    run: ../tools/convert_to_h5.cwl
    in:
      solo_counts_dir: star_solo_align/counts_dir
      raw_count_choice: raw_count_choice
      sample_name: sample_name
      output_basename: output_basename
      soloFeatures: soloFeatures
    out: [raw_converted_h5, filtered_converted_h5]
  tar_solo_count_outdir:
    run: ../tools/tar.cwl
    in:
      output_filename:
        source: output_basename
        valueFrom: $(self).STAR_solo.count.tar.gz
      input_dir: star_solo_align/counts_dir
    out: [output]
  mimic_cr_counts_dir:
    run: ../tools/star_solo_mimic_cr_dir.cwl
    in:
      solo_counts_dir: star_solo_align/counts_dir
      soloFeatures: soloFeatures
      output_name: output_basename
    out: [cr_like_counts_dir]
  tar_solo_cr_mimic_dir:
    run: ../tools/tar.cwl
    in:
      output_filename:
        source: output_basename
        valueFrom: $(self).STAR_Solo_Cell_Ranger_mimic.count.tar.gz
      input_dir: mimic_cr_counts_dir/cr_like_counts_dir
    out: [output]
  seurat_hbc_qc:
    run: ../tools/seurat_hbc_scrna_qc.cwl
    in:
      h5_matrix_inputs: create_h5_output/raw_converted_h5
      sample_id: sample_name
      output_basename: output_basename
      min_umi: qc_min_umi
      min_genes: qc_min_genes
      min_complexity: qc_min_complexity
      max_mito_ratio: qc_max_mito_ratio
      min_gene_prevalence: qc_min_gene_prevalence
    out: [qc_plots, qc_boxplot_stats, cell_counts, seurat_prefilter_data, seurat_filtered_data, variable_features_plot]

sbg:license: Apache License 2.0
sbg:publisher: KFDRC
$namespaces:
  sbg: https://sevenbridges.com
hints:
- class: 'sbg:maxNumberOfParallelInstances'
  value: 2
