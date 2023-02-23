#!/usr/bin/env python

"""
Evaluate integration using cell cycle conservation

Usage:
    metric-cellCycle.py --dataset=<str> --method=<str> --integration=<str> --out-file=<path> <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --out-file=<path>    Path to output file.
"""


def calculate_cellCycleConservation(adata):
    """
    Calculate the cell cycle conservation score for an integrated dataset.

    Parameters
    ----------
    adata
        AnnData object containing the integrated dataset

    Returns
    -------
    Cell cycle conservation score
    """
    from scib.metrics import cell_cycle

    if adata.uns["Species"].lower() not in ["human", "mouse"]:
        from warnings import warn

        warn(
            f"'{adata.uns['Species']}' is not a valid species ('human', 'mouse'). A score of 1 will be returned."
        )
        return 1.0

    print("Calculating cell cycle conservation score...")
    score = cell_cycle(
        adata,
        adata,
        batch_key="Batch",
        embed="X_emb",
        organism=adata.uns["Species"].lower(),
        verbose=True,
    )
    print("Final score: {score}")

    return score


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad
    from functions.functions import format_metric_results

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
    score = calculate_cellCycleConservation(input)
    output = format_metric_results(
        dataset, method, integration, "Integration", "CellCycle", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
