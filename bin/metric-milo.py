#!/usr/bin/env python

"""
Evaluate unseen populations using the MILO metric

Usage:
    metric-milo.py --dataset=<str> --method=<str> --integration=<str> --reference=<path> --out-file=<path> <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --reference=<path>   Path to the reference H5AD file.
    --out-file=<path>    Path to output file.
"""


def calculate_MILO(input):
    """
    Calculate the MILO score for a mapped dataset.

    Parameters
    ----------
    input
        The function input

    Returns
    -------
    The function output
    """

    from scanpy.preprocessing import neighbors
    from milopy.core import make_nhoods, count_nhoods, DA_nhoods

    print("Calculating neighbourhood graph...")
    n_batches = input.obs["Batch"].unique().size
    n_neighbours = min([n_batches * 5, 200])
    neighbors(input, n_neighbors=n_neighbours, use_rep="X_emb")

    print("Creating cell neighbourhoods...")
    # Use neighbourhoods for least 10% or up to 20,000 cells
    prop = max([0.1, min([20000 / input.n_obs, 1])])

    unseen_labels = input.obs["Label"][input.obs["Unseen"]].unique()
    unseen_selected = False
    while not unseen_selected and prop <= 1:
        print(f"Using {prop * input.n_obs:.0f} cells ({prop * 100:.2f}%)")
        make_nhoods(input, prop=prop, seed=1)
        selected_cells = input.obs["nhood_ixs_refined"] == 1
        unseen_counts = (
            input[selected_cells & input.obs["Unseen"]].obs["Label"].value_counts()
        )
        if (unseen_counts.size != unseen_labels.size) or max(unseen_counts) < 10:
            print("Not enough cells in some unseen labels. Selecting more cells...")
            prop = prop * 1.5
        else:
            unseen_selected = True

    print("Counting labels in neighborhoods...")
    count_nhoods(input, sample_col="Batch")

    print("Testing for differential abundance...")
    input.obs["IsQuery"] = input.obs["Dataset"] == "Query"
    DA_nhoods(input, design="~ IsQuery")

    print("Calculating unseen label scores...")
    neighbourhood_adata = input.uns["nhood_adata"]
    label_scores = []
    for label in unseen_labels:
        label_ids = input.obs_names[input.obs["Label"] == label]
        is_label = neighbourhood_adata.obs["index_cell"].isin(label_ids)
        label_fdrs = neighbourhood_adata[is_label].obs["SpatialFDR"]
        # Score is the proportion of neighborhoods with FDR < 0.1
        label_scores.append((label_fdrs < 0.1).mean())

    print("Calculating final MILO score...")
    score = sum(label_scores) / len(label_scores)

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
    score = calculate_MILO(input)
    output = format_metric_results(
        dataset, method, integration, "Unseen", "MILO", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
