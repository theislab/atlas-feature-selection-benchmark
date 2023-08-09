#!/usr/bin/env python

"""
Evaluate integration using Integration LISI (iLISI)

Usage:
    metric-iLISI.py --dataset=<str> --method=<str> --integration=<str> --out-file=<path> <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --out-file=<path>    Path to output file.
"""


def calculate_iLISI(adata):
    """
    Calculate the Integration LISI (iLISI) score for an integrated dataset.

    Parameters
    ----------
    adata
        AnnData object containing the integrated dataset

    Returns
    -------
    The [0, 1] score non-cluster to compact cluster structure.
    """
    from scib.metrics import ilisi_graph
    from scanpy.preprocessing import neighbors
    from igraph import Graph

    print("Calculating final score...")
    score = ilisi_graph(
        adata,
        batch_key="Batch",
        type_="embed",
        use_rep="X_emb",
        k0=90,
        subsample=None,
        scale=True,
        n_cores=2,
        verbose=True,
    )
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
    print("Read data:")
    print(input)
    score = calculate_iLISI(input)
    output = format_metric_results(
        dataset, method, integration, "IntegrationBatch", "iLISI", score
    )
    print(output)
    print("Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
