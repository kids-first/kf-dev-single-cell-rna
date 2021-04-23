cwlVersion: v1.0
class: CommandLineTool
id: hisat2_align
doc: "Run HISAT2 alignment"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0'
  - class: ResourceRequirement
    ramMin: ${return inputs.ram * 1000}
    coresMin: $(inputs.cpus)
  - class: InlineJavascriptRequirement

baseCommand: [tar, xf]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.reference.path) &&

      ${
        var cmd = "hisat2 -t -x " + inputs.reference.nameroot.split('.')[0] +
        "/" + inputs.reference.nameroot.split('.')[0] + " --rg-id=" +
        inputs.input_id + " --rg SM:" + inputs.input_id + " --rg LB:" +
        inputs.input_id + " --rg PL:ILLUMINA --rg PU:" + inputs.input_id +
        " --new-summary --summary-file " + inputs.output_basename +
        ".log --met-file " + inputs.output_basename + ".hisat2.met.txt --met 5 --seed 12345 -k 10 -p " +
        inputs.cpus + " --secondary -S " + inputs.output_basename + "_unsorted.sam "
        if(inputs.fastq2){
          cmd += "-1 " + inputs.fastq1.path + " -2 " + inputs.fastq2.path;
        }
        else {
          cmd += "-U " + inputs.fastq1.path;
        }
        if (inputs.strict){
          cmd += " --mp 1,1 --np 1 --score-min L,0,-0.1 --no-mixed --no-softclip --no-discordant --rdg 99999999,99999999 --rfg 99999999,99999999 --no-spliced-alignment && " +
          "samtools view -S -b " + inputs.output_basename + "_unsorted.sam > " + inputs.output_basename + ".bam"
        }
        else {
          cmd += " && samtools sort -@ " + inputs.cpus + " -O bam -o " +
          inputs.output_basename + ".bam " + inputs.output_basename + "_unsorted.sam && " +
          "samtools index -@ " + inputs.cpus + " " + inputs.output_basename + ".bam " +
          inputs.output_basename + ".bai"
        }
        return cmd;
       }

inputs:
  reference: {type: File, doc: "tarball of reference files"}
  fastq1: {type: File, doc: "gzipped read 1 fq file"}
  fastq2: {type: File?, doc: "gzipped read 2 fq file"}
  output_basename: {type: string, doc: "Output file basename"}
  input_id: {type: string, doc: "Sample name"}
  strict: {type: boolean?, doc: "Flag to use a stricter alignment to make input bam for RSEM"}
  ram: {type: ['null', int], default: 16, doc: "In GB"}
  cpus: {type: ['null', int], default: 4, doc: "Number of CPUs to request"}

outputs:
  log_file:
    type: File
    outputBinding:
      glob: $(inputs.output_basename).log
    doc: "File containing the hisat2 logs"
  met_file:
      type: File
      outputBinding:
        glob: $(inputs.output_basename).hisat2.met.txt
      doc: "File containing the hisat2 metrics"
  bam:
      type: File
      outputBinding:
        glob: $(inputs.output_basename).bam
      doc: "Bam file"
  bai:
      type: File?
      outputBinding:
        glob: $(inputs.output_basename).bai
      doc: "Bam index"
