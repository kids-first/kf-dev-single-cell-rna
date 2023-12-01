cwlVersion: v1.2
class: Workflow
id: kf-single-cell-star-solo-10x-align-wf
label: "KFDRC Single Cell RNA STAR Solo 10x Alignment Workflow"
doc: "STAR Solo Pipeline"
requirements:
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
inputs:
  # multi-step
  output_basename: {type: string, doc: "basename used to name output files"}
  sample_name: {type: 'string', doc: "used as prefix for finding fastqs to analyze,
      e.g. 1k_PBMCs_TotalSeq_B_3p_LT_antibody if the names of the underlying fastqs
      are of the form 1k_PBMCs_TotalSeq_B_3p_LT_antibody_S1_L001_I1_001.fastq.gz,
      one per input fastq in the same order"}
  # star solo
  outSAMattrRGline: { type: string, doc: "Suggested setting, with TABS SEPARATING \
      THE TAGS, format is: ID:sample_name LB:aliquot_id PL:platform SM:BSID for \
      example ID:7316-242 LB:750189 PL:ILLUMINA SM:BS_W72364MN" }
  genomeDir: { type: File, doc: "Tar gzipped reference that will be unzipped at run time" }
  readFilesIn1: { type: 'File[]', doc: "Input fastq file(s), gzipped or uncompressed" }
  readFilesIn2: { type: 'File[]', doc: "R2 or 'mates' reads file(s), gzipped or uncompressed" }
  runThreadN: { type: 'int?', doc: "Num threads for STAR Solo to use", default: 16 }
  solo_type: { type: [ 'null', {type: enum, name: soloType, symbols: ["CB_UMI_Simple", "CB UMI Complex", "CB samTagOut", "SmartSeq"]}],
    doc: "type of single-cell RNAseq. \
    CB_UMI_Simple: (a.k.a. Droplet) one UMI and one Cell Barcode of fixed length in read2, e.g. Drop-seq and 10X Chromium, \
    CB UMI_Complex: multiple Cell Barcodes of varying length, one UMI of fixed length and one adapter sequence of fixed length are allowed in read2 only (e.g. inDrop, ddSeq), \
    CB samTagOut: output Cell Barcode as CR and/or CB SAm tag. No UMI counting. –readFilesIn cDNA read1 [cDNA read2 if paired-end] CellBarcode read . Requires -outSAMtype BAM Unsorted [and/or SortedByCoordinate] \
    Smart-seq: each cell in a separate FASTQ (paired- or single-end), barcodes are corresponding read-groups, no UMI sequences, alignments deduplicated according to alignment start and end (after extending soft-clipped bases)",
    default: "CB_UMI_Simple" }
  soloCBwhitelist: { type: File, doc: "file with whitelist of cell barcodes" }
  soloUMIlen: { type: 'int?', doc: "UMI length", default: 12 }
  clipAdapterType: { type: [ 'null', {type: enum, name: clipAdapterType, symbols: ["Hamming", "CellRanger4", "None"]}], doc: "adapter clipping type. \
    Hamming: adapter clipping based on Hamming distance, with the number of mismatches controlled by -clip5pAdapterMMp. \
    CellRanger4: 5p and 3p adapter clipping similar to CellRanger4. Utilizes Opalpackage by Martin Sosic: https://github.com/Martinsos/opal,\
    None: no adapter clipping, all other clip* parameters are disregarded",
    default: "CellRanger4" }
  outFilterScoreMin: { type: 'int?', doc: "alignment will be output only if its score is higher than or equal to this value", default: 30 }
  soloCBmatchWLtype: { type: [ 'null', {type: enum, name: clipAdapterType, symbols: ["Exact", "1MM", "1MM_multi", "1MM_multi_pseudocounts", "1MM_multi_Nbase_pseudocounts", "EditDist_2"]}], doc: "matching the Cell Barcodes to the WhiteList. \
    Exact: only exact matches allowed. \
    1MM: only one match in whitelist with 1 mismatched base allowed. Allowed CBs have to have at least one read with exact match. \
    1MM_multi: multiple matches in whitelist with 1 mismatched base allowed, posterior probability calculation is used choose one of the matches. Allowed CBs have to have at least one read with exact match. This option matches best with CellRanger 2.2.0. \
    1MM_multi_pseudocounts: same as 1MM Multi, but pseudocounts of 1 are added to all whitelist barcodes. \
    1MM_multi_Nbase_pseudocounts: same as 1MM multi pseudocounts, multimatching to WL is allowed for CBs with N-bases. This option matches best with CellRanger >= 3.0.0 \
    EditDist_2: allow up to edit distance of 3 fpr each of the barcodes. May include one deletion + one insertion. Only works with -soloType CB UMI Complex. Matches to multiple passlist barcdoes are not allowed. Similar to ParseBio Split-seq pipeline.",
    default: "1MM_multi_Nbase_pseudocounts" }
  soloUMIfiltering: { type: [ 'null', {type: enum, name: soloUMIfiltering, symbols: ["-", "MultiGeneUMI", "MultiGeneUMI_All", "MultiGeneUMI_CR"]}], doc: "type of UMI filtering (for reads uniquely mapping to genes). -: basic filtering: remove UMIs with N and homopolymers (similar to CellRanger 2.2.0). \
    MultiGeneUMI: basic + remove lower-count UMIs that map to more than one gene. \
    MultiGeneUMI_All: basic + remove all UMIs that map to more than one gene. \
    MultiGeneUMI_CR: basic + remove lower-count UMIs that map to more than one gene, matching CellRanger > 3.0.0. Only works with -soloUMIdedup 1MM_CR",
    default: "MultiGeneUMI_CR" }
  soloUMIdedup:  { type: [ 'null', {type: enum, name: soloUMIdedup, symbols: ["1MM_All", "1MM_Directional_UMItools", "1MM_Directional", "Exact", "NoDedup", "1MM_CR"]}], doc: "type of UMI deduplication (collapsing) algorithm). \
    1MM_All: all UMIs with 1 mismatch distance to each other are collapsed (i.e. counted once). \
    1MM_Directional_UMItools: follows the ”directional” method from the UMI-tools by Smith, Heger and Sudbery (Genome Research 2017). \
    1MM_Directional: same as 1MM Directional UMItools, but with more stringent criteria for duplicate UMIs. \
    Exact: only exactly matching UMIs are collapsed \
    NoDedup: no deduplication of UMIs, count all reads. \
    1MM_CR: CellRanger2-4 algorithm for 1MM UMI collapsing.",
    default: "1MM_CR" }
  soloFeatures: { type: [ 'null', {type: enum, name: soloFeatures, symbols: ["Gene", "SJ", "GeneFull", "GeneFull_ExonOverIntron", "GeneFull_Ex50pAS"]}], doc: "genomic features for which the UMI counts per Cell Barcode are collected. \
    Gene: genes: reads match the gene transcript. \
    SJ: splice junctions: reported in SJ.out.tab. \
    GeneFull: full gene (pre-mRNA): full gene (pre-mRNA): count all reads overlapping genes' exons and introns. \
    GeneFull_ExonOverIntron: count all reads overlapping genes' exons and introns: prioritize 100% overlap with exons. \
    GeneFull_Ex50pAS: full gene (pre-RNA): count all reads overlapping genes' exons and introns: prioritize >50% overlap with exons. Do not count reads with 100% exonic overlap in the antisense direction.",
    default: "GeneFull" }
  soloCellFilter: { type: [ 'null', {type: enum, name: soloCellFilter, symbols: ["None", "TopCells", "CellRanger2.2", "EmptyDrops_CR"]}], doc: "cell filtering type and parameters. \
    None: do not output filtered cells. \
    TopCells: only report top cells by UMI count, followed by the exact number of cells. \
    CellRanger2.2: simple filtering of CellRanger 2.2. Can be followed by numbers: number of expected cells, robust maximum percentile for UMI count, maximum to minimum ratio for UMI count. The harcoded values are from CellRanger: nExpectedCells=3000; maxPercentile=0.99; maxMinRatio=10. \
    EmptyDrops_CR: EmptyDrops filtering in CellRanger flavor. Please cite the original EmptyDrops paper: A.T.L Lun et al, Genome Biology, 20, 63 (2019): https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y, Can be followed by 10 numeric parameters: nExpectedCells maxPercentile maxMinRatio indMin indMax umiMin umiMinFracMedian candMaxN FDR simN. The harcoded values are from CellRanger: 3000 0.99 10 45000 90000 500 0.01 20000 0.01 10000",
    default: "EmptyDrops_CR" }
  outSAMtype: { type: [ 'null', {type: enum, name: outSAMtype, symbols: ["BAM Unsorted", "None", "BAM SortedByCoordinate", "SAM Unsorted", "SAM SortedByCoordinate"]}],
    default: "BAM SortedByCoordinate",
    doc: "type of SAM/BAM output. None: no SAM/BAM output. Otherwise, first word is output type (BAM or SAM), second is sort type (Unsorted or SortedByCoordinate)" }
  # seurat qc
  seurat_qc_min_genes: {type: "int?", doc: "minimum number of genes per cell", default: 400}
  seurat_qc_max_genes: {type: "int?", doc: "maximum number of genes per cell", default: 4000}
  seurat_qc_max_mt: {type: "float?", doc: "maximum percent mitochondrial reads per cell",
    default: 5}
  seurat_qc_normalize_method: {type: ['null', {type: enum, name: normalize_method,
        symbols: ["log_norm", "sct"]}], default: "log_norm", doc: "normalization method.
      One of log_norm or sct"}
  seurat_qc_num_pcs: {type: "int?", doc: "number of PCs to calculate", default: 30}
outputs:
  star_solo_counts_dir: {type: File, outputSource: tar_solo_count_outdir/output }
  star_solo_bam: {type: 'File?', outputSource: star_solo_align/genomic_bam_out }
  star_solo_matrix_filtered: {type: File, outputSource: create_filtered_h5_output/converted_h5 }
  star_solo_matrix_raw: {type: File, outputSource: create_raw_h5_output/converted_h5 }
  star_solo_log_final_out: { type: File, outputSource: star_solo_align/log_final_out }
  star_solo_junctions: { type: File, outputSource: star_solo_align/junctions_out }
  seurat_qc_html: {type: File, outputSource: rename_seurat_html/renamed_file }
  seurat_qc_rds: {type: File, outputSource: rename_seurat_rds/renamed_file }
steps:
  star_solo_align:
    run: ../tools/star_solo_2.7.10b.cwl
    in:
      outSAMattrRGline: outSAMattrRGline
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
      outSAMtype: outSAMtype
    out: [log_progress_out, log_out, log_final_out, genomic_bam_out, junctions_out, counts_dir]
  create_raw_h5_output:
    run: ../tools/convert_to_h5.cwl
    in:
      # dir structure is parent/Gene*/raw
      counts_matrix_dir:
        source: star_solo_align/counts_dir
        valueFrom: |
          ${for(var i in self.listing){
            var item = self.listing[i];
            if (item.basename.startsWith("Gene")){
            for(var j in item.listing){
                var subItem = item.listing[j];
                if (subItem.basename == "raw"){
                    return subItem;
                    }
                }
              }
            }
          }
      sample_name: sample_name
      output_basename: output_basename
    out: [converted_h5]
  create_filtered_h5_output:
    run: ../tools/convert_to_h5.cwl
    in:
      # dir structure is parent/Gene*/filtered
      counts_matrix_dir:
        source: star_solo_align/counts_dir
        valueFrom: |
          ${for(var i in self.listing){
            var item = self.listing[i];
            if (item.basename.startsWith("Gene")){
            for(var j in item.listing){
                var subItem = item.listing[j];
                if (subItem.basename == "filtered"){
                    return subItem;
                    }
                }
              }
            }
          }
      sample_name: sample_name
      output_basename: output_basename
    out: [converted_h5]
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
      solo_counts_dir:
        source: star_solo_align/counts_dir
        valueFrom: |
          ${for(var i in self.listing){
            var item = self.listing[i];
            if (item.basename.startsWith("Gene")){
              return item;
              }
            }
          }
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
  seurat_qc:
    run: ../tools/seurat_qc.cwl
    in:
      filtered_bc_matrix_dir: mimic_cr_counts_dir/cr_like_counts_dir
      sample_name: sample_name
      min_genes: seurat_qc_min_genes
      max_genes: seurat_qc_max_genes
      max_mt: seurat_qc_max_mt
      normalize_method: seurat_qc_normalize_method
      num_pcs: seurat_qc_num_pcs
    out: [result_dir, summary_html, rds]
  rename_seurat_html:
    run: ../tools/rename_file.cwl
    in:
      in_file: seurat_qc/summary_html
      out_filename:
        source: output_basename
        valueFrom: $(self).seurat.qc.html
    out: [renamed_file]
  rename_seurat_rds:
    run: ../tools/rename_file.cwl
    in:
      in_file: seurat_qc/rds
      out_filename:
        source: output_basename
        valueFrom: $(self).seurat.qc.rds
    out: [renamed_file]
sbg:license: Apache License 2.0
sbg:publisher: KFDRC
$namespaces:
  sbg: https://sevenbridges.com
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 2
