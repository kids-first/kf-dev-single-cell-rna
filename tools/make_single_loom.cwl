cwlVersion: v1.0
class: CommandLineTool
id: make_single_loom
doc: "Create a single sample loom file from the rsem gene output."

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/loompy:2.0.16'
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: make_rsem_loom.py
        entry: |-
          '''Script to create count and QC looms.'''
          import sys
          import os
          import csv
          import argparse
          import loompy
          import numpy as np
          import scipy as sc
          def parse_args(args):
              '''Get arguments.'''
              parser = argparse.ArgumentParser()
              #required args
              required_args = parser.add_argument_group("required arguments")
              required_args.add_argument(
                  "--data",
                  help="Path to the RSEM output file.",
                  required=True
                  )
              required_args.add_argument(
                  "--sample",
                  help="Sample name.",
                  required=True
                  )
              args = parser.parse_args()
              data_file = args.data
              sample_name = args.sample
              if not os.path.isfile(data_file):
                  raise ValueError(data_file + " is not a file or does not exist.")
              return data_file, sample_name
          def make_count_matrix(counts):
              '''Convert counts list to matrix.'''
              nrows, ncols = np.shape(counts)
              coo = sc.sparse.coo_matrix(counts[:])
              x = []
              y = []
              value = []
              for k in range(0, coo.data.shape[0]):
                  x.append(coo.row[k])
                  y.append(coo.col[k])
                  value.append(coo.data[k])
              x = np.asarray(x)
              y = np.asarray(y)
              value = np.asarray(value)
              matrix = sc.sparse.coo_matrix((value, (y, x)), shape=(ncols, nrows))
              return matrix
          def make_loom(file, sample_name):
              '''Turn the RSEM output into loom file, Separate gene id and names.'''
              reader = csv.DictReader(open(file), delimiter="\t")
              count_metric = 'TPM'
              ens_ids = []
              genes = []
              counts = []
              for line in reader:
                  gene = line['gene_id']
                  if "_" in gene:
                      gene_id, name = gene.split("_",1)
                  counts.append(float(line[count_metric]))
                  ens_ids.append(gene_id)
                  genes.append(name)
              count_mat = make_count_matrix([counts])
              loomfile = sample_name + ".loom"
              col_attrs = dict()
              col_attrs['sample_name'] = np.array([sample_name])
              row_attrs = {"ensemble_ids":np.array(ens_ids), "gene_names":np.array(genes)}
              loompy.create(loomfile, count_mat, row_attrs, col_attrs)
          def main(args):
              '''Main, read data file and make loom.'''
              data_file, sample_name = parse_args(args)
              make_loom(data_file, sample_name)
          if __name__ == "__main__":
              # execute only if run as a script
              main(sys.argv)

        writable: false

baseCommand: [python]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      make_rsem_loom.py --data $(inputs.rsem_matrix.path) --sample $(inputs.sample_name)

inputs:
  rsem_matrix: {type: 'File', doc: "Gene level output file from rsem"}
  sample_name: {type: string, doc: "Will be used as the sample's name inside the loom file and as the basename for the output loom file."}

outputs:
  loom_file:
    type: File
    outputBinding:
      glob: $(inputs.sample_name).loom
    doc: "Output loom file."
