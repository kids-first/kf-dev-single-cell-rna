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
  help="Comma separated list of paths to the RSEM output looms",
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
