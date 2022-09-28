#!/usr/bin/env python

"""
Evaluate integration using isolated label scores

Usage:
    metric-isolatedLabels.py --dataset=<str> --method=<str> --integration=<str> --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --cluster            If set use the clustering-based F1 score, otherwise use silhouette scores.
    --out-file=<path>    Path to output file.
"""


def calculate_isolatedLabels(adata, cluster):
    """
    Calculate the Isolated label score for an integrated dataset.

    Parameters
    ----------
    adata
        AnnData object containing the integrated dataset
    cluster
        Whether to use the clustering-based F1 score (True), rather than the silhouette score (False)

    Returns
    -------
    The mean isolated label score
    """
    from scib.metrics import isolated_labels

    print(f"Calculating isolated labels score (cluster={cluster})...")
    score = isolated_labels(adata, 'Label', 'Batch', 'X_emb', cluster=cluster, verbose=True)
    print("Final score: {score}")

    return score


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad
    from _functions import format_metric_results

    args = docopt(__doc__)

    file = args["<file>"]
    dataset = args["--dataset"]
    method = args["--method"]
    integration = args["--integration"]
    cluster = args["--cluster"]
    out_file = args["--out-file"]

    print("Reading data from '{file}'...")
    input = read_h5ad(file)
    print("Read data:")
    print(input)
    score = calculate_isolatedLabels(input, cluster)
    if cluster:
        metric = "IsolatedLabelF1"
    else:
        metric = "IsolatedLabelASW"
    output = format_metric_results(
        dataset, method, integration, "Integration", metric, score
    )
    print(output)
    print("Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
