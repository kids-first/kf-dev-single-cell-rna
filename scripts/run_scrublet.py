'''Compare a new study's QC metrics to a set of studies metrics.'''
import argparse
import sys
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scrublet as scr

def parse_args(args):
    '''Get arguments.'''
    #optional args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-e", "--expected_doublet_rate",
        help='Optional, expected doublet rate, usually specific to the method; \
        default 0.06 for 10X',
        type = float,
        default = 0.06,
        required=False
        )
    parser.add_argument(
        "-s", "--score",
        help='Optional, doublet cut-off, cells with greater scores will be \
        labelled as doublets; must be between 0 and 1; default 0.25',
        type = float,
        default = 0.25,
        required=False
        )
    parser.add_argument(
        "-c", "--min_count",
        help='Optional, minimum expression count to retain a gene; default 2',
        type = int,
        default = 2,
        required=False
        )
    parser.add_argument(
        "-l", "--min_cell",
        help='Optional, minimum number of cells a gene must be in to be \
        reatained; default 3',
        type = int,
        default = 3,
        required=False
        )
    parser.add_argument(
        "-g", "--min_gene_variability_pctl",
        help='Optional, Keep the most highly variable genes \
        (in the top min_gene_variability_pctl percentile), as measured by \
        the v-statistic; default 85',
        type = int,
        default = 85,
        required=False
        )
    parser.add_argument(
        "-p", "--n_prin_comps",
        help='Optional, Number of PCs to use for clustering; default 30',
        type = int,
        default = 30,
        required=False
        )

    #required args
    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument(
        "--output",
        help="Output basename",
        required=True
        )
    required_args.add_argument(
        "--matrix",
        help="Input matrix file",
        required=True
        )

    #parse and return arguments
    args = parser.parse_args()
    matrix = args.matrix
    output = args.output
    edr = args.expected_doublet_rate
    score = args.score
    counts = args.min_count
    cells = args.min_cell
    variability = args.min_gene_variability_pctl
    pcs = args.n_prin_comps

    return matrix, output, edr, score, counts, cells, variability, pcs

def main(args):
    '''Main, take args, run script.'''
    mat_dir, out_base, edr, score_thresh, counts, cells, variability, pcs \
        = parse_args(args)

    #read matrix file
    count_mat = sc.read_10x_mtx(mat_dir, cache=False)

    #score doublets
    scrub = scr.Scrublet(count_mat, expected_doublet_rate=edr)
    count_mat.obs['doublet_scores'], count_mat.obs['predicted_doublets'] = \
        scrub.scrub_doublets(min_counts = counts,
        min_cells = cells, min_gene_variability_pctl = variability,
        n_prin_comps = pcs, use_approx_neighbors = False)
    count_mat.obs['predicted_doublets'] = scrub.call_doublets(threshold=score_thresh)

    #plot scores histogram
    hist_file = out_base + ".hist.png"
    scrub.plot_histogram()
    plt.savefig(hist_file)

    #write output file with doublet scores and predicted doublets
    score_out_file = out_base + ".csv"
    pd.DataFrame(count_mat.obs).to_csv(score_out_file)

if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv)
