cwlVersion: v1.0
class: ExpressionTool
id: build_samples_array
doc: "Turn a file with sample names into an array of strings."
requirements:
  - class: InlineJavascriptRequirement

inputs:
  sample_file:
    type: 'File'
    doc: "File with sample names"
    inputBinding:
      loadContents: true

outputs:
  sample_names:
    type: string[]

expression:
  "${
    var unfiltered = inputs.sample_file.contents.split('\\n');
    var sample_names = unfiltered.filter(function(el) { return el; });
    return {'sample_names': sample_names}
  }"
