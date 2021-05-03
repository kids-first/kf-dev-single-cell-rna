cwlVersion: v1.0
class: CommandLineTool
id: merge_looms
doc: "Merge several loom files into one"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/loompy:2.0.16'
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: assemble_ss2_output.py
        entry: |-
          '''Script to assemble count and QC data from individual looms.'''
          import sys
          import re
          import os
          import argparse
          import loompy

          def parse_args(args):
              '''Get arguments.'''
              #required args
              parser = argparse.ArgumentParser()
              required_args = parser.add_argument_group("required arguments")
              required_args.add_argument(
                  "--files",
                  help="Comma separated list of aths to the directory with the RSEM output looms",
                  required=True
                  )
              required_args.add_argument(
                  "--base",
                  help="Output basename",
                  required=True
                  )

              #parse and return arguments
              args = parser.parse_args()
              files = args.files
              base = args.base

              return files, base

          def main(args):
              '''Main.'''
              #parse arguments
              files, base_name = parse_args(args)
              merged_loom = base_name + ".loom"

              #make matrix loom
              loom_files = files.split(',')
              print(loom_files)
              loompy.combine(loom_files, output_file=merged_loom)

          if __name__ == "__main__":
              # execute only if run as a script
              main(sys.argv)

        writable: false

baseCommand: []

arguments:
  - position: 1
    shellQuote: false
    valueFrom: |-
     locs="${
       var arr = [];
       for (var i=0; i<inputs.loom_files.length; i++)
        arr = arr.concat(inputs.loom_files[i].path)
       return (arr.join(' '))
     }"
     looms="${
       var arr = [];
       for (var i=0; i<inputs.loom_files.length; i++)
        arr = arr.concat("/tmp/" + inputs.loom_files[i].basename)
       return (arr.join(','))
     }"
     cp -t /tmp/ $locs
     python assemble_ss2_output.py --base $(inputs.output_basename) --files $looms

inputs:
  output_basename: {type: string, doc: "basename of output file"}
  loom_files: {type: 'File[]', doc: "List of loom files to combine"}

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_basename).loom
    doc: "Combined matrix loom file"
