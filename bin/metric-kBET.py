#!/usr/bin/env python

"""
Evaluate integration using the k-Nearest Neighbour Batch effect Test (k-BET) score.
The k-BET score measures the mixing (distribution) of batches and labels within the neighbourhood of a cell.
This is to evaluate how well represented each batch is in a given label.

Usage:
    metric-kBET.py --dataset=<str> --method=<str> --integration=<str> --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --out-file=<path>    Path to output file.
"""


def calculate_kBET(adata):
    """
    Calculate the k-Nearest Neighbour Batch effect Test (k-BET) score for an
    integrated dataset over batches per group(e. g., cell type).

    Parameters
    ----------
    adata
        AnnData object containing the integrated dataset

    Returns
    -------
    The scaled average of the rejection rates per label; 0 means optimal batch mixing and 1 means low batch mixing.
    """
    from scib.metrics import kBET

    print("Calculating score...")
    score = kBET(
        adata,
        batch_key="Batch",
        label_key="Label",
        type_="embed",
        embed="X_emb",
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
    score = calculate_kBET(input)
    output = format_metric_results(
        dataset, method, integration, "Integration", "kBET", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
