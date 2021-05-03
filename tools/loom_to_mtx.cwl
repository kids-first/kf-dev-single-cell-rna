cwlVersion: v1.0
class: CommandLineTool
id: loom_to_mtx
doc: "Convert a loom file to a Matrix Market file that Seurat can read."

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/loompy:2.0.16'
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: loom_to_mtx.py
        entry: |-
          #!/usr/bin/env python
          import sys
          import h5py
          import numpy as np
          import scipy.io
          import scipy.sparse

          loom_file = sys.argv[1]
          out_dir = sys.argv[2]

          loom = h5py.File(loom_file)

          sample_names = np.array(loom['col_attrs']['sample_name'])
          ensembl_names = np.array(loom['row_attrs']['ensemble_ids'])
          hugo_names = np.array(loom['row_attrs']['gene_names'])
          matrix = np.array(loom['matrix'])

          loom.close()

          #output files
          out_mat = out_dir + '/matrix.mtx'
          out_bar = out_dir + '/barcodes.tsv'
          out_feat = out_dir + '/features.tsv'

          print(out_dir)
          print(out_mat)

          scipy.io.mmwrite(out_mat, scipy.sparse.csr_matrix(matrix), field='real', precision=2)

          with open(out_bar, 'w') as barcodes:
              for name in sample_names[:-1]:
                  barcodes.write(name.decode() + '\n')
              barcodes.write(sample_names[-1].decode())

          with open(out_feat, 'w') as features:
              for ensbl, hugo in list(zip(ensembl_names, hugo_names))[:-1]:
                  features.write('\t'.join([ensbl.decode(), hugo.decode(), 'Gene Expression']) + '\n')
              features.write('\t'.join([ensembl_names[-1].decode(), hugo_names[-1].decode(), 'Gene Expression']))


        writable: false

baseCommand: [mkdir, -p]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      $(inputs.output_basename) &&
      python loom_to_mtx.py $(inputs.loom_file.path) $(inputs.output_basename) &&
      gzip $(inputs.output_basename)/* &&
      tar cvzf $(inputs.output_basename).tar.gz $(inputs.output_basename)

inputs:
  loom_file: {type: 'File', doc: "Loom file to convert"}
  output_basename: {type: string, doc: "Basename for the output tarball"}

outputs:
  output_tar:
    type: File
    outputBinding:
      glob: $(inputs.output_basename).tar.gz
    doc: "Output tarball"
