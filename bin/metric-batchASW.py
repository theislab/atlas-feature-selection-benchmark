#!/usr/bin/env python

"""
Evaluate integration using a Modified Adjusted Silhouette Width (mASW).
The silhouette score is calculated over batches per a group, e. g., cell type.
This is to evaluate how well represented each batch is in a given group.

Usage:
    metric-batchASW.py --dataset=<str> --method=<str> --integration=<str> --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --out-file=<path>    Path to output file.
"""


def calculate_batch_asw(adata):
    """
    Calculate a modified Adjusted Silhouette Width (mASW) score for an integrated dataset over batches per group(e. g., cell type).

    Parameters
    ----------
    adata
        AnnData object containing the integrated dataset

    Returns
    -------
    The scaled mASW score [0, 1] where 1 indicates optimal batch mixing on average across groups (e. g., cell type).
    """
    from scib.metrics import silhouette_batch

    print("Calculating final score...")
    score = silhouette_batch(adata, batch_key="Batch", label_key="Label", embed="X_emb")
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
    score = calculate_batch_asw(input)
    output = format_metric_results(
        dataset, method, integration, "Integration", "batchASW", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
