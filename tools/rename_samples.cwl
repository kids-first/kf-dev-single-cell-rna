cwlVersion: v1.2
class: CommandLineTool
id: cellranger_count
doc: "Run cellranger count on a set of fastq files"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/cellranger:6.0'
  - class: ResourceRequirement
    ramMin: 20000
    coresMin: 16
  - class: InlineJavascriptRequirement
  - class: LoadListingRequirement
    loadListing: shallow_listing

baseCommand: [mkdir, -p]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
     $(inputs.sample_name) &&
     ${
      var cp_cmd = "";
      var files = inputs.fastqs.listing;
      for (var i = 0; i < files.length; i++) {
        var file = files[i];
        if (file.basename.includes(inputs.sample_name)) {
          if (file.basename.includes("R1") || file.basename.includes("r1")) {
            if (inputs.corrected_read_1_name) {
              cp_cmd = cp_cmd + "cp " + inputs.fastqs.path + "/" + file.basename + " " + inputs.sample_name + "/" + inputs.corrected_read_1_name + ".fastq.gz && ";
            }
            else {
              cp_cmd = cp_cmd + "cp " + inputs.fastqs.path + "/" + file.basename + " " + inputs.sample_name + "/" + file.basename + " && ";
            }
          }
          if (file.basename.includes("R2") || file.basename.includes("r2")) {
            if (inputs.corrected_read_2_name) {
              cp_cmd = cp_cmd + "cp " + inputs.fastqs.path + "/" + file.basename + " " + inputs.sample_name + "/" + inputs.corrected_read_2_name + ".fastq.gz && ";
            }
            else {
              cp_cmd = cp_cmd + "cp " + inputs.fastqs.path + "/" + file.basename + " " + inputs.sample_name + "/" + file.basename + " && ";
            }
          }
        }
      }
      return cp_cmd;
     }
     echo "Done"

inputs:
  fastqs: {type: Directory, doc: "directory of fastqs being run"}
  sample_name: {type: string, doc: "sample name, used as prefix for finding fastqs to analyze"}
  corrected_read_1_name: {type: "string?", doc: "corrected read one name"}
  corrected_read_2_name: {type: "string?", doc: "corrected read two name"}

outputs:
  renamed_dir:
    type: Directory
    outputBinding:
      glob: $(inputs.sample_name)
    doc: "Directory containing renamed fastq files"