#!/usr/bin/env python

"""
Evaluate integration using principal component regression

Usage:
    metric-batchPCR.py --dataset=<str> --method=<str> --integration=<str> --out-file=<path> <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --out-file=<path>    Path to output file.
"""


def calculate_batchPCR(adata):
    """
    Calculate the principal component regression score for an integrated dataset.

    Parameters
    ----------
    adata
        AnnData object containing the integrated dataset

    Returns
    -------
    The PCR comparison score
    """
    from scib.metrics import pcr_comparison

    print("Calculating final score...")
    score = pcr_comparison(adata, adata, "Batch", embed="X_emb", n_comps=50, scale=True, verbose=True)
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
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    print("Read data:")
    print(input)
    score = calculate_batchPCR(input)
    output = format_metric_results(
        dataset, method, integration, "Integration", "BatchPCR", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
