cwlVersion: v1.2
class: CommandLineTool
id: star_solo_2-7-10b_align
label: "STAR Solo 2.7.10b Alignment Tool"
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/star:2.7.10b'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.runThreadN)
    ramMin: 60000

baseCommand: [tar, -I pigz, -xvf]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.genomeDir.path)
  - position: 2
    shellQuote: false
    valueFrom: >-
      && STAR
      --genomeDir ./$(inputs.genomeDir.nameroot.replace(".tar", ""))/
      --readFilesCommand $(inputs.readFilesIn1[0].nameext == '.gz' ? 'zcat' : '-')
      --outFileNamePrefix $(inputs.outFileNamePrefix).

inputs:
  genomeDir: { type: File, doc: "Tar gzipped reference that will be unzipped at run time" }
  readFilesIn1: { type: 'File[]', doc: "Input fastq file(s), gzipped or uncompressed",
    inputBinding: { itemSeparator: ",", separate: true, position: 4} }
  readFilesIn2: { type: 'File[]', doc: "R2 or 'mates' reads file(s), gzipped or uncompressed",
    inputBinding: { prefix: "--readFilesIn", itemSeparator: ",", position: 3} }
  runThreadN: { type: 'int?', default: 16, inputBinding: { position: 5, prefix: '--runThreadN' } }
  outFileNamePrefix: { type: string, doc: "output files name prefix (including full or relative path). Can only be defined on the command line. \
    Tool will add '.' after prefix to easily delineate between file name and suffix" }
  twopassMode: { type: ['null', {type: enum, name: twopassMode, symbols: ["Basic", "None"]}], default: "None",
    doc: "Enable two pass mode to detect novel splice events. Default is None (off).", inputBinding: { position: 5, prefix: '--twopassMode' } }
  outSAMattrRGline: { type: 'string?', doc: "Set if outputting bam, with TABS SEPARATING \
      THE TAGS, format is: ID:sample_name LB:aliquot_id PL:platform SM:BSID for \
      example ID:7316-242 LB:750189 PL:ILLUMINA SM:BS_W72364MN",
      inputBinding: { position: 5, prefix: '--outSAMattrRGline', shellQuote: false } }
  outSAMattributes: { type: 'string?', doc: "a string of desired SAM attributes, in the order desired for the output SAM. Tags can be listed in any combination/order. \
      Please refer to the STAR manual, as there are numerous combinations: https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf",
      inputBinding: { position: 5, prefix: '--outSAMattributes', shellQuote: false } }
  solo_type: { type: [ 'null', {type: enum, name: soloType, symbols: ["CB_UMI_Simple", "CB UMI Complex", "CB samTagOut", "SmartSeq"]}],
    doc: "type of single-cell RNAseq. \
    CB_UMI_Simple: (a.k.a. Droplet) one UMI and one Cell Barcode of fixed length in read2, e.g. Drop-seq and 10X Chromium, \
    CB UMI_Complex: multiple Cell Barcodes of varying length, one UMI of fixed length and one adapter sequence of fixed length are allowed in read2 only (e.g. inDrop, ddSeq), \
    CB samTagOut: output Cell Barcode as CR and/or CB SAm tag. No UMI counting. –readFilesIn cDNA read1 [cDNA read2 if paired-end] CellBarcode read . Requires -outSAMtype BAM Unsorted [and/or SortedByCoordinate] \
    Smart-seq: each cell in a separate FASTQ (paired- or single-end), barcodes are corresponding read-groups, no UMI sequences, alignments deduplicated according to alignment start and end (after extending soft-clipped bases)",
    default: "CB_UMI_Simple",
    inputBinding: { position: 5, prefix: "--soloType", shellQuote: false } }
  soloCBwhitelist: { type: File, doc: "file with whitelist of cell barcodes",
     inputBinding: { position: 5, prefix: "--soloCBwhitelist"} }
  soloUMIlen: { type: 'int?', doc: "UMI length", default: 12,
    inputBinding: { position: 5, prefix: "--soloUMIlen"} }
  clipAdapterType: { type: [ 'null', {type: enum, name: clipAdapterType, symbols: ["Hamming", "CellRanger4", "None"]}], doc: "adapter clipping type. \
    Hamming: adapter clipping based on Hamming distance, with the number of mismatches controlled by -clip5pAdapterMMp. \
    CellRanger4: 5p and 3p adapter clipping similar to CellRanger4. Utilizes Opalpackage by Martin Sosic: https://github.com/Martinsos/opal,\
    None: no adapter clipping, all other clip* parameters are disregarded",
    default: "CellRanger4",
    inputBinding: { position: 5, prefix: "--clipAdapterType"} }
  outFilterScoreMin: { type: 'int?', doc: "alignment will be output only if its score is higher than or equal to this value", default: 30,
    inputBinding: { position: 5, prefix: "--outFilterScoreMin"} }
  soloCBmatchWLtype: { type: [ 'null', {type: enum, name: clipAdapterType, symbols: ["Exact", "1MM", "1MM_multi", "1MM_multi_pseudocounts", "1MM_multi_Nbase_pseudocounts", "EditDist_2"]}], doc: "matching the Cell Barcodes to the WhiteList. \
    Exact: only exact matches allowed. \
    1MM: only one match in whitelist with 1 mismatched base allowed. Allowed CBs have to have at least one read with exact match. \
    1MM_multi: multiple matches in whitelist with 1 mismatched base allowed, posterior probability calculation is used choose one of the matches. Allowed CBs have to have at least one read with exact match. This option matches best with CellRanger 2.2.0. \
    1MM_multi_pseudocounts: same as 1MM Multi, but pseudocounts of 1 are added to all whitelist barcodes. \
    1MM_multi_Nbase_pseudocounts: same as 1MM multi pseudocounts, multimatching to WL is allowed for CBs with N-bases. This option matches best with CellRanger >= 3.0.0 \
    EditDist_2: allow up to edit distance of 3 fpr each of the barcodes. May include one deletion + one insertion. Only works with -soloType CB UMI Complex. Matches to multiple passlist barcdoes are not allowed. Similar to ParseBio Split-seq pipeline.",
    default: "1MM_multi_Nbase_pseudocounts",
    inputBinding: { position: 5, prefix: "--soloCBmatchWLtype"} }
  soloUMIfiltering: { type: [ 'null', {type: enum, name: soloUMIfiltering, symbols: ["-", "MultiGeneUMI", "MultiGeneUMI_All", "MultiGeneUMI_CR"]}], doc: "type of UMI filtering (for reads uniquely mapping to genes). -: basic filtering: remove UMIs with N and homopolymers (similar to CellRanger 2.2.0). \
    MultiGeneUMI: basic + remove lower-count UMIs that map to more than one gene. \
    MultiGeneUMI_All: basic + remove all UMIs that map to more than one gene. \
    MultiGeneUMI_CR: basic + remove lower-count UMIs that map to more than one gene, matching CellRanger > 3.0.0. Only works with -soloUMIdedup 1MM_CR",
    default: "MultiGeneUMI_CR",
    inputBinding: { position: 5, prefix: "--soloUMIfiltering"} }
  soloUMIdedup:  { type: [ 'null', {type: enum, name: soloUMIdedup, symbols: ["1MM_All", "1MM_Directional_UMItools", "1MM_Directional", "Exact", "NoDedup", "1MM_CR"]}], doc: "type of UMI deduplication (collapsing) algorithm). \
    1MM_All: all UMIs with 1 mismatch distance to each other are collapsed (i.e. counted once). \
    1MM_Directional_UMItools: follows the ”directional” method from the UMI-tools by Smith, Heger and Sudbery (Genome Research 2017). \
    1MM_Directional: same as 1MM Directional UMItools, but with more stringent criteria for duplicate UMIs. \
    Exact: only exactly matching UMIs are collapsed \
    NoDedup: no deduplication of UMIs, count all reads. \
    1MM_CR: CellRanger2-4 algorithm for 1MM UMI collapsing.",
    default: "1MM_CR",
    inputBinding: { position: 5, prefix: "--soloUMIdedup"} }
  soloFeatures: { type: [ 'null', {type: enum, name: soloFeatures, symbols: ["Gene", "SJ", "GeneFull", "GeneFull_ExonOverIntron", "GeneFull_Ex50pAS"]}], doc: "genomic features for which the UMI counts per Cell Barcode are collected. \
    Gene: genes: reads match the gene transcript. \
    SJ: splice junctions: reported in SJ.out.tab. \
    GeneFull: full gene (pre-mRNA): full gene (pre-mRNA): count all reads overlapping genes' exons and introns. \
    GeneFull_ExonOverIntron: count all reads overlapping genes' exons and introns: prioritize 100% overlap with exons. \
    GeneFull_Ex50pAS: full gene (pre-RNA): count all reads overlapping genes' exons and introns: prioritize >50% overlap with exons. Do not count reads with 100% exonic overlap in the antisense direction.",
    default: "GeneFull",
    inputBinding: { position: 5, prefix: "--soloFeatures"} }
  soloCellFilter: { type: [ 'null', {type: enum, name: soloCellFilter, symbols: ["None", "TopCells", "CellRanger2.2", "EmptyDrops_CR"]}], doc: "cell filtering type and parameters. \
    None: do not output filtered cells. \
    TopCells: only report top cells by UMI count, followed by the exact number of cells. \
    CellRanger2.2: simple filtering of CellRanger 2.2. Can be followed by numbers: number of expected cells, robust maximum percentile for UMI count, maximum to minimum ratio for UMI count. The harcoded values are from CellRanger: nExpectedCells=3000; maxPercentile=0.99; maxMinRatio=10. \
    EmptyDrops_CR: EmptyDrops filtering in CellRanger flavor. Please cite the original EmptyDrops paper: A.T.L Lun et al, Genome Biology, 20, 63 (2019): https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y, Can be followed by 10 numeric parameters: nExpectedCells maxPercentile maxMinRatio indMin indMax umiMin umiMinFracMedian candMaxN FDR simN. The harcoded values are from CellRanger: 3000 0.99 10 45000 90000 500 0.01 20000 0.01 10000",
    default: "EmptyDrops_CR",
    inputBinding: { position: 5, prefix: "--soloCellFilter"} }
  soloMultiMappers: { type: 'string[]?', doc: "Possible one or more values: 'Unique', 'Uniform', 'PropUnique', 'EM', 'Rescue'. Including multi-gene reads allows for more accurate gene quantification and, more importantly, enables detection of gene expression from certain classes of genes that are supported only by multi-gene reads, such as overlapping genes and highly similar paralog families. \
    Unique: software default, count only reads that map to unique genes \
    Uniform: uniformly distributes the multi-gene UMIs to all genes in its gene set. Each gene gets a fractional count of 1/N_genes, where N_genes is the number of genes in the set. This is the simplest possible option, and it offers higher sensitivity for gene detection at the expense of lower precision \
    PropUnique: distributes the multi-gene UMIs proportionally to the number of unique UMIs per gene. UMIs that map to genes that are not supported by unique UMIs are distributed uniformly \
    EM: uses Maximum Likelihood Estimation (MLE) to distribute multi-gene UMIs among their genes, taking into account other UMIs (both unique- and multi-gene) from the same cell (i.e. with the same CB). Expectation-Maximization (EM) algorithm is used to find the gene expression values that maximize the likelihood function. Recovering multi-gene reads via MLE-EM model was previously used to quantify transposable elements in bulk RNA-seq {TEtranscripts} and in scRNA-seq {Alevin; Kallisto-bustools}. \
    Rescue: distributes multi-gene UMIs to their gene set proportionally to the sum of the number of unique-gene UMIs and uniformly distributed multi-gene UMIs in each gene Mortazavi et al. It can be thought of as the first step of the EM algorithm",
    inputBinding: { position: 5, prefix: "--soloMultiMappers" } }
  extra_args: { type: 'string?', inputBinding: { position: 5, shellQuote: false }, doc: "Any additional arguments for this tool. See STAR Documentation for complete list of options. Example input: --limitSjdbInsertNsj 1000001" }
  outSAMtype: { type: [ 'null', {type: enum, name: outSAMtype, symbols: ["BAM Unsorted", "None", "BAM SortedByCoordinate", "SAM"]}],
    default: "None",
    doc: "type of SAM/BAM output. None: no SAM/BAM output. Otherwise, first word is output type (BAM or SAM), second, if BAM, is sort type (Unsorted or SortedByCoordinate)",
    inputBinding: { position: 3, prefix: '--outSAMtype', shellQuote: false } }

outputs:
  log_progress_out: { type: File, doc: "Simple progress output. Can use to gauge speed and run time", outputBinding: {glob: '*Log.progress.out'} }
  log_out: { type: File, doc: "Contains a summary of all params used and reference files", outputBinding: {glob: '*Log.out'} }
  log_final_out: { type: File, doc: "Overall summary of read mapping statistics", outputBinding: {glob: '*Log.final.out'} }
  genomic_bam_out: { type: 'File?', doc: "Reads mapping to genomic coordinates", outputBinding: {glob: '*Aligned*[b|s]am'} }
  junctions_out: { type: File, doc: "high confidence collapsed splice junctions in tab-delimited form", outputBinding: {glob: '*SJ.out.tab'} }
  counts_dir: { type: Directory, doc: "Output dir with raw and filtered counts", outputBinding: { glob: '*.Solo.out'} }
