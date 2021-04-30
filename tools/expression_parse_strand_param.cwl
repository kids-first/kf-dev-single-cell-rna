cwlVersion: v1.0
class: ExpressionTool
id: expression_strand_params
requirements:
  - class: InlineJavascriptRequirement

inputs:
  wf_strand_param: {type: [{type: enum, name: wf_strand_param, symbols: ["default", "rf-stranded", "fr-stranded"]}], doc: "use 'default' for unstranded/auto, rf_stranded if read1 in the fastq read pairs is reverse complement to the transcript, fr-stranded if read1 same sense as transcript"}

outputs:
  rsem_std:
    type: string
  rnaseqc_std:
    type: string

expression:
  "${
      var strand = 'default';
      if (inputs.wf_strand_param != null){
        strand = inputs.wf_strand_param;
      }
      var parse_dict = {
          'default': {'rsem_std': 'none', 'rnaseqc_std': 'default'},
          'rf-stranded': {'rsem_std': 'reverse', 'rnaseqc_std': 'rf'},
          'fr-stranded': {'rsem_std': 'forward', 'rnaseqc_std': 'fr'}
          };
      if (strand in parse_dict){
        return parse_dict[strand];

      }
      else{
        throw new Error(strand + ' is a not a valid strand param. Use one of default, rf-stranded, fr-stranded');
      }
  }"
