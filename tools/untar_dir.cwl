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
outputs:
  outdir: 
    type: Directory
    outputBinding:
      glob: $(inputs.tarfile.nameroot.split('.')[0])
      loadListing: deep_listing
