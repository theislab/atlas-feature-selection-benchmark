#!/usr/bin/env python

"""
Evaluate mapping using Query LISI (qLISI)

Usage:
    metric-qLISI.py --dataset=<str> --method=<str> --integration=<str> --out-file=<path> <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --out-file=<path>    Path to output file.
"""


def calculate_qLISI(adata):
    """
    Calculate the Query LISI (qLISI) score for a mapped dataset.

    Parameters
    ----------
    adata
        AnnData object containing the mapped dataset

    Returns
    -------
    The [0, 1] score for query batch mixing.
    """
    from scib.metrics import ilisi_graph
    from scanpy.preprocessing import neighbors
    from igraph import Graph

    # Return 1 if there is only one batch
    if adata.obs["Batch"].nunique() == 1:
        import warnings

        warnings.warn("Only one batch, returning a score of 1")
        return 1

    print("Calculating qLISI score...")
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
    score = calculate_qLISI(input)
    output = format_metric_results(
        dataset, method, integration, "Mapping", "qLISI", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
