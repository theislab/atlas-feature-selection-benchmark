#!/usr/bin/env python

"""
Evaluate integration using graph connectivity

Usage:
    metric-graphConnectivity.py --dataset=<str> --method=<str> --integration=<str> --out-file=<path> <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --out-file=<path>    Path to output file.
"""


def calculate_graphConnectivity(adata):
    """
    Calculate the graph connectivity score for an integrated dataset.

    Parameters
    ----------
    adata
        AnnData object containing the integrated dataset

    Returns
    -------
    Graph connectivity score
    """
    from scib.metrics import graph_connectivity
    from scanpy.preprocessing import neighbors

    print("Calculating nearest neighbors...")
    neighbors(adata, use_rep="X_emb")

    print("Calculating graph connectivity score...")
    score = graph_connectivity(adata, label_key="Label")
    print("Final score: {score}")

    return score


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad
    from functions.metrics import format_metric_results

    args = docopt(__doc__)

    file = args["<file>"]
    dataset = args["--dataset"]
    method = args["--method"]
    integration = args["--integration"]
    out_file = args["--out-file"]

    print("Reading data from '{file}'...")
    input = read_h5ad(file)
    input.obs["Label"] = input.obs["Label"].cat.remove_unused_categories()
    print("Read data:")
    print(input)
    score = calculate_graphConnectivity(input)
    output = format_metric_results(
        dataset, method, integration, "IntegrationBatch", "GraphConnectivity", score
    )
    print(output)
    print("Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
