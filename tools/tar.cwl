cwlVersion: v1.2
class: CommandLineTool
id: tar
requirements:
  - class: InlineJavascriptRequirement
  - class: LoadListingRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.cpu)
    ramMin: $(inputs.ram * 1000)
  - class: InitialWorkDirRequirement
    listing: $(inputs.input_dir)
baseCommand: [tar, czf]
inputs:
  output_filename:
    type: string
    inputBinding:
      position: 1
  input_dir:
    type: Directory
    loadListing: deep_listing
    inputBinding:
      position: 2
      valueFrom: $(self.basename)
  # Control
  cpu: { type: 'int?', default: 4, doc: "Number of threads to use." }
  ram: { type: 'int?', default: 16, doc: "GB of RAM to allocate to this task." }
outputs:
  output: 
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
