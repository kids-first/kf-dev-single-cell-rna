cwlVersion: v1.2
class: CommandLineTool
id: rnaseqc
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'gcr.io/broad-cga-aarong-gtex/rnaseqc:2.4.2'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cpus)
    ramMin: ${return inputs.ram * 1000}

baseCommand: [rnaseqc]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.collapsed_gtf.path)
      $(inputs.input_bam.path)
      output/
      ${
        var cmd = "";
        if (inputs.legacy) {
          cmd += " --legacy"
        }
        if (inputs.strand != null && inputs.strand != "default"){
          cmd += " --stranded=" + inputs.strand;
        }
        if (! inputs.paired) {
          cmd += " --unpaired";
        }
        return cmd;
      }

inputs:
  input_bam: {type: File, secondaryFiles: [^.bai]}
  collapsed_gtf: {type: File, doc: "GTF file without overlapping genes and one transcript per gene"}
  strand: {type: 'string?', doc: "Only collect metrics for features on the same strand. Input string must be: RF, rf, FR, or fr"}
  paired: {type: 'boolean?', default: TRUE}
  legacy: {type: 'boolean?', default: FALSE, doc: "use legacy RNAseQC counting rules. Legacy is default counting scheme used by KF bulk RNA workflow"}
  ram: {type: 'int?', default: 8, doc: "In GB"}
  cpus: {type: 'int?', default: 4, doc: "Number of CPUs to request"}

outputs:
  metrics:
    type: File
    outputBinding:
      glob: 'output/*.metrics.tsv'
  gene_TPM:
    type: File
    outputBinding:
      glob: 'output/*.gene_tpm.gct'
  gene_count:
    type: File
    outputBinding:
      glob: 'output/*.gene_reads.gct'
  exon_count:
    type: File
    outputBinding:
      glob: 'output/*.exon_reads.gct'
