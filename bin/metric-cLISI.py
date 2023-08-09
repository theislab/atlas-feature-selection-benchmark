#!/usr/bin/env python

"""
Evaluate integration using Cell-type LISI (cLISI)

Usage:
    metric-cLISI.py --dataset=<str> --method=<str> --integration=<str> --out-file=<path> <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --out-file=<path>    Path to output file.
"""


def calculate_cLISI(adata):
    """
    Calculate the Cell-type LISI (cLISI) score for an integrated dataset.

    Parameters
    ----------
    adata
        AnnData object containing the integrated dataset

    Returns
    -------
    cLISI score
    """

    from scib.metrics import clisi_graph
    from scanpy.preprocessing import neighbors
    from igraph import Graph

    adata.obs["Label"] = adata.obs["Label"].cat.remove_unused_categories()

    print("Calculating cLISI score...")
    score = clisi_graph(
        adata,
        label_key="Label",
        type_="embed",
        use_rep="X_emb",
        k0=90,
        subsample=None,
        scale=True,
        n_cores=2,
        verbose=True,
    )
    print(f"Final score: {score}")

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

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    print("Read data:")
    print(input)
    score = calculate_cLISI(input)
    output = format_metric_results(
        dataset, method, integration, "IntegrationBio", "cLISI", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
