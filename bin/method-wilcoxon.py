#!/usr/bin/env python

"""
Select Wilcoxon marker genes

Usage:
    method-wilcoxon.py --out-file=<path> [options] <file>

Options:
    -h --help               Show this screen.
    -o --out-file=<path>    Path to output file.
    -n --n-features=<int>   Number of features to select per label [default: 200].
"""


def select_features_wilcoxon(adata, n_features):
    """
    Select features using the Wilcoxon rank sum test to detect marker genes

    Parameters
    ----------
    adata
        AnnData object
    n_features
        Number of features to select per label

    Returns
    ----------
    DataFrame containing the selected features
    """

    from scanpy.preprocessing import normalize_total, log1p
    from scanpy.tools import rank_genes_groups, filter_rank_genes_groups
    from scanpy.get import rank_genes_groups_df
    from pandas import DataFrame

    print("Normalising expression...")
    normalize_total(adata, target_sum=1e4)
    log1p(adata)

    print("Testing for marker genes...")
    # Remove labels with no cells (unseen labels)
    adata.obs["Label"] = adata.obs["Label"].cat.remove_unused_categories()
    rank_genes_groups(adata, groupby="Label", method="wilcoxon", tie_correct=True)

    print("Filtering marker genes...")
    filter_rank_genes_groups(
        adata, min_in_group_fraction=0.1, max_out_group_fraction=0.8
    )

    selected_features = []
    for label in adata.obs["Label"].cat.categories:
        print(f"Getting results for label '{label}'...")
        filtered_results = rank_genes_groups_df(
            adata, group=str(label), key="rank_genes_groups_filtered"
        )
        filtered_results = filtered_results[filtered_results["names"].notnull()]

        filtered_results = filtered_results[filtered_results["pvals_adj"] <= 0.01]

        filtered_results = filtered_results.sort_values(
            by="logfoldchanges", ascending=False
        )
        filtered_results = filtered_results.head(n=n_features)

        selected_features = selected_features + list(filtered_results["names"])

    selected_features = list(set(selected_features))
    print(f"Selected {len(selected_features)} features...")
    output = DataFrame(selected_features, columns=["Feature"])

    return output


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad

    args = docopt(__doc__)

    file = args["<file>"]
    n_features = int(args["--n-features"])
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    print("Read data:")
    print(input)
    output = select_features_wilcoxon(input, n_features)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
