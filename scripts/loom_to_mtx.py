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

scipy.io.mmwrite(out_mat, scipy.sparse.csr_matrix(matrix), field='real', precision=2)

with open(out_bar, 'w') as barcodes:
    for name in sample_names:
        barcodes.write(name.decode() + '\n')

with open(out_feat, 'w') as features:
    for ensbl, hugo in list(zip(ensembl_names, hugo_names)):
        features.write('\t'.join([ensbl.decode(), hugo.decode(), 'Gene Expression']) + '\n')
