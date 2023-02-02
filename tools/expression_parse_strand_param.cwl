cwlVersion: v1.2
class: ExpressionTool
id: expression_strand_params
doc: "Format the overall strand parameter to the format expected by downstream tools (rsem and RNAseQC)"
requirements:
  - class: InlineJavascriptRequirement

inputs:
  wf_strand_param: {type: [{type: enum, name: wf_strand_param, symbols: ["default", "rf-stranded", "fr-stranded"]}], doc: "use 'default' for unstranded/auto, rf_stranded if read1 in the fastq read pairs is reverse complement to the transcript, fr-stranded if read1 same sense as transcript"}

outputs:
  rsem_std:
    type: string
  rnaseqc_std:
    type: string
  hisat2_std: { type: ['null', {type: enum, name: rna_strandness, symbols: ["FR", "RF", "F", "R"]}] }

expression:
  "${
      var parse_dict = {
          'default': {'rsem_std': 'none', 'rnaseqc_std': 'default'},
          'rf-stranded': {'rsem_std': 'reverse', 'rnaseqc_std': 'rf', 'hisat2_std': 'RF'},
          'fr-stranded': {'rsem_std': 'forward', 'rnaseqc_std': 'fr', 'hisat2_std': 'FR'}
          };
      return parse_dict[inputs.wf_strand_param];
  }"
