#!/usr/bin/env python

"""
Evaluate the cluster structure of the integrated data using the Normalized Mutual Information (NMI).
This is a similarity metrix between "ground truth" clusters and found clusters after integration.
This score is not adjusted for chance.

Usage:
    metric-nmi.py --dataset=<str> --method=<str> --integration=<str> --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --out-file=<path>    Path to output file.
"""


def calculate_nmi(adata):
    """
    Calculate the Normalized Mutual Information (NMI) score for a integrated dataset.

    Parameters
    ----------
    adata
        AnnData object containing the integrated dataset

    Returns
    -------
    The NMI score [0, 1] where 0 is no mutual information, and 1 a perfect correlation.
    """
    from scib.metrics import nmi
    from scib.metrics import opt_louvain

    print("Optimising clusters...")
    opt_louvain(adata, label_key="Label", cluster_key="Cluster", function=nmi)
    print("Calculating score...")
    score = nmi(adata, group1="Label", group2="Cluster")
    print(f"Final score: {score}")

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
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    print("Read data:")
    print(input)
    score = calculate_nmi(input)
    output = format_metric_results(
        dataset, method, integration, "Integration", "NMI", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()