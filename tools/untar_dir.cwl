cwlVersion: v1.2
class: CommandLineTool
id: untar_dir
requirements:
  - class: InitialWorkDirRequirement
    listing: $(inputs.tarfile)
baseCommand: [tar, -zxf]
inputs:
  tarfile:
    type: File
    inputBinding:
      position: 1
  untar_dir_name: { type: 'string?', doc: "If untar dir name different than expected file nameroot, provide here" }
outputs:
  outdir: 
    type: Directory
    outputBinding:
      glob: |
        $(inputs.untar_dir_name == null ? inputs.tarfile.nameroot.split('.')[0] : inputs.untar_dir_name)
      loadListing: deep_listing
