#!/usr/bin/env python

"""
Evaluate integration using principal component regression

Usage:
    metric-batchPCR.py --dataset=<str> --method=<str> --integration=<str> --exprs=<file> --out-file=<path> <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --exprs=<file>       Path to H5AD file containing the expression matrix.
    --out-file=<path>    Path to output file.
"""


def calculate_batchPCR(adata, exprs):
    """
    Calculate the principal component regression score for an integrated dataset.

    Parameters
    ----------
    adata
        AnnData object containing the integrated dataset
    exprs
        AnnData object containing the unintegrated dataset with the expression matrix

    Returns
    -------
    The PCR comparison score
    """
    from scib.metrics import pcr_comparison

    print("Calculating final score...")
    score = pcr_comparison(
        adata_pre=exprs,
        adata_post=adata,
        covariate="Batch",
        embed="X_emb",
        n_comps=50,
        scale=True,
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
    exprs_file = args["--exprs"]
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    print("Read data:")
    print(input)
    print(f"Reading expression data from '{exprs_file}'...")
    exprs = read_h5ad(exprs_file)
    print("Read expression data:")
    print(exprs)
    score = calculate_batchPCR(input, exprs)
    output = format_metric_results(
        dataset, method, integration, "IntegrationBatch", "BatchPCR", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
