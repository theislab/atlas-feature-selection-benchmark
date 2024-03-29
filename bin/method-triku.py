#!/usr/bin/env python

"""
Select features using the Triku method

Usage:
    method-triku.py --out-file=<path> <file>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
"""


def select_triku_features(adata):
    """
    Select features in a dataset using Triku

    Parameters
    ----------
    adata
        AnnData object

    Returns
    -------
    DataFrame containing the selected features
    """

    import triku as tk
    import scanpy as sc

    print("Filtering and normalizing data...")
    sc.pp.filter_cells(adata, min_genes=50)
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    print("Building neighbourhood graph...")
    sc.pp.pca(adata)
    sc.pp.neighbors(adata, metric="cosine", n_neighbors=int(0.5 * len(adata) ** 0.5))

    print("Selecting triku features...")
    tk.tl.triku(adata, verbose="info")

    adata.var["Feature"] = adata.var.index
    selected_features = adata.var
    selected_features = selected_features[selected_features["highly_variable"] == True]

    return selected_features


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad

    args = docopt(__doc__)

    file = args["<file>"]
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    print("Read data:")
    print(input)
    output = select_triku_features(input)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
