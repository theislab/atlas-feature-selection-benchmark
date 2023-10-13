#!/usr/bin/env python

"""
Evaluate integration using Query-Reference Mapping LISI (mLISI)

Usage:
    metric-mLISI.py --dataset=<str> --method=<str> --integration=<str> --reference=<path> --out-file=<path> <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --reference=<path>   Path to the reference H5AD file.
    --out-file=<path>    Path to output file.
"""


def calculate_mLISI(adata):
    """
    Calculate the Mapping LISI (mLISI) score for a mapped dataset.

    Parameters
    ----------
    adata
        AnnData object containing the mapped dataset

    Returns
    -------
    The [0, 1] score for query-reference mixing.
    """
    from scib.metrics import ilisi_graph

    print("Calculating mLISI score...")
    # We use the iLISI function but with query/reference as the batch key
    score = ilisi_graph(
        adata,
        batch_key="Dataset",
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

    query_file = args["<file>"]
    dataset = args["--dataset"]
    method = args["--method"]
    integration = args["--integration"]
    reference_file = args["--reference"]
    out_file = args["--out-file"]

    print(f"Reading query data from '{query_file}'...")
    query = read_h5ad(query_file)
    query.obs["Dataset"] = "Query"
    print("Read data:")
    print(query)
    print(f"Reading reference data from '{reference_file}'...")
    reference = read_h5ad(reference_file)
    reference.obs["Dataset"] = "Reference"
    print("Read data:")
    print(reference)
    print("Merging query and reference...")
    input = reference.concatenate(query)
    input.obs["Dataset"] = input.obs["Dataset"].astype("category")
    print("Merged data:")
    print(input)
    score = calculate_mLISI(input)
    output = format_metric_results(
        dataset, method, integration, "Mapping", "mLISI", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
