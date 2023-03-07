#!/usr/bin/env python

"""
Select features using the hotspot method

Usage:
    method-hotspot.py --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
    --n-features=<int>   Number of features to select [default: 1000].
"""


def select_hotspot_features(adata, n):
    """
    Select HotSpot features

    Parameters
    ----------
    adata
        AnnData object
    n
        Number of features to select

    Returns
    -------
    DataFrame containing the selected features
    """

    from pandas import DataFrame
    from scipy.sparse import csc_matrix
    import scanpy as sc
    import hotspot

    adata.obs["total_counts"] = adata.X.sum(axis=1)
    # Save counts for use by HotSpot
    # Tutorial recommends a CSC matrix for improved speed
    adata.layers["counts"] = csc_matrix(adata.X.copy())

    print("Calculating PCA...")
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.scale(adata)
    sc.tl.pca(adata)
    sc.pp.pca(adata)

    print("Calculating HotSpot features...")
    hs = hotspot.Hotspot(
        adata,
        layer_key="counts",
        model="danb",
        latent_obsm_key="X_pca",
        umi_counts_obs_key="total_counts",
    )
    hs.create_knn_graph(weighted_graph=False, n_neighbors=30)
    hs_results = hs.compute_autocorrelations()

    print("Selecting top features...")
    # Filtering like done in tutorial (not really for feature selection, more for specific downstream analyses):
    # https://hotspot.readthedocs.io/en/latest/CD4_Tutorial.html#Grouping-genes-into-lineage-based-modules
    selected_features = (
        hs_results.loc[hs_results["FDR"] < 0.05]
        .sort_values("Z", ascending=False)
        .head(n)
        .index
    )

    selected_features = DataFrame({"Feature": selected_features})

    return selected_features


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
    output = select_hotspot_features(input, n_features)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
