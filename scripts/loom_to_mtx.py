#!/usr/bin/env python

import sys

import h5py
import numpy as np
import scipy.io
import scipy.sparse

matrix_header = '%%MatrixMarket matrix coordinate integer general\n%metadata_json: {"format_version": 2, "software_version": "3.1.0"}\n'

loom_file = sys.argv[1]

loom = h5py.File(loom_file)

sample_names = np.array(loom['col_attrs']['sample_name'])
ensembl_names = np.array(loom['row_attrs']['ensemble_ids'])
hugo_names = np.array(loom['row_attrs']['gene_names'])
matrix = np.array(loom['matrix'])

loom.close()

total_cells = len(sample_names)
total_features = len(hugo_names)
total_counts = np.sum(matrix)

scipy.io.mmwrite('matrix.mtx', scipy.sparse.csr_matrix(matrix), field='real', precision=2)

with open('barcodes.tsv', 'w') as barcodes:
    for name in sample_names[:-1]:
        barcodes.write(name.decode() + '\n')
    barcodes.write(sample_names[-1].decode())

with open('features.tsv', 'w') as features:
    for ensbl, hugo in list(zip(ensembl_names, hugo_names))[:-1]:
        features.write('\t'.join([ensbl.decode(), hugo.decode(), 'Gene Expression']) + '\n')
    features.write('\t'.join([ensembl_names[-1].decode(), hugo_names[-1].decode(), 'Gene Expression']))
