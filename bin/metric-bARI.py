#!/usr/bin/env python

"""
Evaluate the cluster structure of the integrated data using the balanced Adjusted Rand Index (bARI).
This is a similarity metrix between "ground truth" clusters and found clusters after integration adjusted for class imbalance.

Usage:
    metric-bARI.py --dataset=<str> --method=<str> --integration=<str> --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --out-file=<path>    Path to output file.
"""


def calculate_bARI(adata):
    """
    Calculate the balanced Adjusted Rand Index (bARI) score for a integrated dataset.

    Parameters
    ----------
    adata
        AnnData object containing the integrated dataset

    Returns
    -------
    The bARI score [0, 1] where 0 is no match between the pair of labels, and 1 a perfect match.
    """
    from scib.metrics import cluster_optimal_resolution
    from balanced_clustering import balanced_adjusted_rand_index
    from scanpy.preprocessing import neighbors

    print("Calculating neighbors...")
    neighbors(adata, use_rep="X_emb")
    print("Optimising clusters...")
    cluster_optimal_resolution(
        adata, label_key="Label", cluster_key="Cluster", metric=bARI
    )
    print("Calculating score...")
    score = balanced_adjusted_rand_index(
        labels_true=adata.obs["Label"], labels_pred=adata.obs["Cluster"]
    )
    print(f"Final score: {score}")

    return score


def bARI(adata, label_key, cluster_key):
    """
    Balanced Adjusted Rand Index (bARI) for an AnnData object.
    A wrapper around `balanced_clustering.balanced_adjusted_rand_index()` with the interface expected by `scib.metrics.cluster_optimal_resolution`.

    Parameters
    ----------
    adata
        Anndata object with cluster assignments
    cluster_key
        Name of the column in adata.obs containing cluster assignments
    label_key
        Name of the column in adata.obs containing ground truth labels

    Returns
    -------
    The bARI score [0, 1] where 0 is no match between the pair of labels, and 1 a perfect match.
    """

    from balanced_clustering import balanced_adjusted_rand_index

    clusters = adata.obs[cluster_key].to_numpy()
    labels = adata.obs[label_key].to_numpy()

    return balanced_adjusted_rand_index(labels_true=labels, labels_pred=clusters)


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad
    from functions.functions import format_metric_results

    args = docopt(__doc__)

    file = args["<file>"]
    dataset = args["--dataset"]
    method = args["--method"]
    integration = args["--integration"]
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    print("Read data:")
    print(input)
    score = calculate_bARI(input)
    output = format_metric_results(
        dataset, method, integration, "Integration", "bARI", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
