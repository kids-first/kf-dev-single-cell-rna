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

    #parse and return arguments
    args = parser.parse_args()
    data_file = args.data
    sample_name = args.sample

    if not os.path.isfile(data_file):
        raise ValueError(data_file + " is not a file or does not exist.")

    return data_file, sample_name

def make_count_matrix(counts):
    '''Convert counts list to matrix.'''
    #this might flips the matrix from one row with many columns to many colums, one row

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

    #define which RSEM output field to use as count; assuming this should be TPM
    #but I might be wrong and / or we may want this to be configurable in future
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

    #convert counts to matrix
    count_mat = make_count_matrix([counts])

    #make loom
    loomfile = sample_name + ".loom"
    col_attrs = dict()
    col_attrs['sample_name'] = np.array([sample_name])
    row_attrs = {"ensemble_ids":np.array(ens_ids), "gene_names":np.array(genes)}
    loompy.create(loomfile, count_mat, row_attrs, col_attrs)

def main(args):
    '''Main, read data file and make loom.'''

    #parse arguments
    data_file, sample_name = parse_args(args)

    #make matrix loom
    make_loom(data_file, sample_name)

if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv)
