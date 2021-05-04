cwlVersion: v1.0
class: ExpressionTool
id: build_fastq2_array
doc: "If the fastq2 parameter isn't given, build an array of empty strings that are the same length as the fastq1s input."
requirements:
  - class: InlineJavascriptRequirement

inputs:
  fastq1s: {type: 'File[]', doc: "Array of fastq 1s to align"}
  fastq2s: {type: 'File[]?', doc: "Array of fastq 2s to align"}

outputs:
  fastq2_array:
    type: File[]

expression:
  "${
      var fastq2s = [];
      if (inputs.fastq2s == null){
        for (var i=0; i<inputs.fastq1s.length; i++)
          fastq2s = fastq2s.concat(null)
      }
      else {
        fastq2s = inputs.fastq2s;
      }
    return {'fastq2_array': fastq2s}
  }"
