cwlVersion: v1.2
class: CommandLineTool
id: cellranger_mkref
doc: "Make a custom reference for cellranger. Adjusting the attributes field will subset the gtf file to limit the reference to include only genes meeting that criteria"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/cellranger:6.1.2'
  - class: ResourceRequirement
    ramMin: $(inputs.cr_instance_ram * 1000)
    coresMin: $(inputs.nthreads)
  - class: InlineJavascriptRequirement

baseCommand: []

arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
     ${
        var cmd = "";
        if (inputs.attribute != null){
            cmd = "cellranger mkgtf " + inputs.genes.path + " " + inputs.genes.nameroot + ".filtered.gtf";
        }
        return cmd;
     }
  - position: 1
    shellQuote: false
    valueFrom: >-
     ${
        var cmd = "";
        if (inputs.attribute != null){
            cmd = " && cellranger mkref --genes=" + inputs.genes.nameroot + ".filtered.gtf";
        }
        else{
            cmd = "cellranger mkref --genes=" + inputs.genes.path;
        }
        return cmd;
     }

inputs:
  nthreads: { type: 'int?', doc: "Num cores to use", default: 16,
    inputBinding: {position: 2, prefix: "--nthreads=", separate: false } }
  cr_instance_ram: { type: 'int?', doc: 'Ram in GB to make available to STAR index step', default: 64,
    inputBinding: { position: 2, prefix: "--memgb=", separate: false} }
  genome: { type: string, doc: "genome name", 
    inputBinding: { position: 2, prefix: "--genome=", separate: false } }
  fasta: { type: File, doc: "Reference fasta file to index", 
    inputBinding: { position: 2, prefix: "--fasta=", separate: false } }
  genes: { type: File, doc: "GTF file with gene model to use" }
  attribute: { type: ["null", { type: array, items: string, inputBinding: {prefix: "--attribute=", separate: false }}], doc: "allowed biotypes", 
    inputBinding: { position: 0, separate: false } }

outputs:
  cellranger_reference:
    type: Directory
    outputBinding:
      glob: $(inputs.genome)
    doc: "Cell ranger ref dir"
