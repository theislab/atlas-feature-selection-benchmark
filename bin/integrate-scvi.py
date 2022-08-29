#!/usr/bin/env python

"""
Integrate a dataset using scVI

Usage:
    integrate-scvi.py --features=<path> --out-dir=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --features=<path>    Path to a file containing selected features.
    --out-dir=<path>     Path to output directory.
"""

import scvi
from _functions import add_umap, add_integrated_embeddings, suffix_embeddings, plot_embedding


def run_scVI(adata):
    """
    Integrate a dataset using scVI

    Parameters
    ----------
    adata
        AnnData object containing the dataset to integrate

    Returns
    -------
    Integrated scVI model
    """

    print("Setting up AnnData for scVI...")
    scvi.model.SCVI.setup_anndata(adata, batch_key="Batch")
    print(adata)

    print("Creating scVI model...")
    # Recommended parameters when you plan to use scArches reference mapping
    # from https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scarches_scvi_tools.html
    arches_params = dict(
        use_layer_norm="both",
        use_batch_norm="none",
        encode_covariates=True,
        dropout_rate=0.2,
        n_layers=2,
    )
    model = scvi.model.SCVI(adata, **arches_params)
    print(model)
    model.view_anndata_setup()

    print("Training scVI model...")
    model.train()
    print(model)

    return model


def select_features(adata, features):
    """
    Subset a dataset to selected features

    Parameters
    ----------
    adata
        AnnData object containing the dataset to subset
    features
        DataFrame containing the selected features

    Returns
    -------
    AnnData containing only the selected features
    """

    print(f"Subsetting to {features.shape[0]} selected features...")
    adata = adata[:, features["Feature"]].copy()

    return adata


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad
    from pandas import read_csv
    from os.path import join

    args = docopt(__doc__)

    file = args["<file>"]
    features_file = args["--features"]
    out_dir = args["--out-dir"]

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    print("Read data:")
    print(input)
    print(f"Reading selected features from '{features_file}'...")
    features = read_csv(features_file, sep="\t")
    print("Read features:")
    print(features)
    input = select_features(input, features)
    output = run_scVI(input)
    print("Adding unintegrated UMAP...")
    add_umap(output.adata)
    suffix_embeddings(output.adata)
    add_integrated_embeddings(output, output.adata)
    print(f"Writing output to '{out_dir}'...")
    output.save(out_dir, save_anndata=True, overwrite=True)
    umap_file = join(out_dir, "umap-unintegrated.png")
    umap = plot_embedding(output.adata, basis="X_umap_unintegrated")
    print(f"Saving unintegrated UMAP plot to '{umap_file}'...")
    umap.savefig(umap_file, dpi=300, bbox_inches="tight")
    umap_file = join(out_dir, "umap-integrated.png")
    umap = plot_embedding(output.adata)
    print(f"Saving integrated UMAP plot to '{umap_file}'...")
    umap.savefig(umap_file, dpi=300, bbox_inches="tight")
    print("Done!")


if __name__ == "__main__":
    main()
