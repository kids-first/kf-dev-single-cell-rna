cwlVersion: v1.2
class: CommandLineTool
id: rsem-calculate-expression
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'images.sbgenomics.com/uros_sipetic/rsem:1.3.1'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cpus)
    ramMin: ${return inputs.ram * 1000}

baseCommand: [tar]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      -zxf $(inputs.reference.path) &&
      ${
        var cmd = "rsem-calculate-expression --alignments --append-names --no-bam-output -p " +
        inputs.cpus + " --bam --seed 12345 --calc-pme --single-cell-prior"
        if (inputs.paired) {
          cmd += " --paired-end";
        }
        cmd += " " + inputs.input_bam.path +
        " ./" + inputs.reference.nameroot.split('.')[0] + "/" + inputs.reference.nameroot.split('.')[0] +
        " " +  inputs.output_basename;
        return cmd
      }

inputs:
  input_bam: {type: File}
  reference: {type: File, doc: "tarball of reference files"}
  output_basename: {type: string, doc: "Output file basename"}
  paired: {type: ['null', boolean], default: TRUE}
  ram: {type: ['null', int], default: 24, doc: "In GB"}
  cpus: {type: ['null', int], default: 16, doc: "Number of CPUs to request"}

outputs:
  gene_out:
    type: File
    outputBinding:
      glob: $(inputs.output_basename).genes*

  isoform_out:
    type: File
    outputBinding:
      glob: $(inputs.output_basename).isoforms*

  cnt_out:
    type: File
    outputBinding:
      glob: $(inputs.output_basename).stat/$(inputs.output_basename).cnt

  model_out:
    type: File
    outputBinding:
      glob: $(inputs.output_basename).stat/$(inputs.output_basename).model

  theta_out:
    type: File
    outputBinding:
      glob: $(inputs.output_basename).stat/$(inputs.output_basename).theta
