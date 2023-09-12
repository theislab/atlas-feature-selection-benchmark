#!/usr/bin/env python

"""
Integrate a dataset using Symphony

Usage:
    integrate-symphony.py --features=<path> --out-dir=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --features=<path>    Path to a file containing selected features.
    --out-dir=<path>     Path to output directory.
    --seed=<int>         Random seed to use [default: 0].
"""

from functions.integration import (
    add_umap,
    suffix_embeddings,
    plot_embedding,
)


def run_symphony(adata, seed):
    """
    Integrate a dataset using Symphony

    Parameters
    ----------
    adata
        AnnData object containing the dataset to integrate
    seed
        Random seed to use

    Returns
    -------
    AnnData object with integrated embedding in adata.obsm["X_emb"]
    """

    from scanpy.preprocessing import scale
    from scanpy.tools import pca
    from symphonypy.preprocessing import harmony_integrate

    print("Scaling...")
    scale(adata, max_value=10)

    print("Calculating PCA...")
    pca(adata, n_comps=30, zero_center=False)

    print("Integrating with Harmony...")
    print(f"Using random seed: {seed}")
    harmony_integrate(
        adata, key="Batch", ref_basis_adjusted="X_emb", verbose=True, random_seed=seed
    )

    return adata


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad
    from scanpy.preprocessing import normalize_total, log1p
    from pandas import read_csv
    from os.path import join
    from functions.anndata import minimise_anndata
    import os

    args = docopt(__doc__)

    file = args["<file>"]
    features_file = args["--features"]
    out_dir = args["--out-dir"]
    seed = int(args["--seed"])

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    print("Read data:")
    print(input)

    print("Storing counts matrix...")
    counts = input.X.copy()

    print("Normalising...")
    normalize_total(input, target_sum=1e4)

    print("Log transforming...")
    log1p(input)

    print(f"Reading selected features from '{features_file}'...")
    features = read_csv(features_file, sep="\t")
    print("Read features:")
    print(features)

    print(f"Subsetting to {features.shape[0]} selected features...")
    input = input[:, input.var_names.isin(features["Feature"])].copy()
    counts = counts[:, input.var_names.isin(features["Feature"])].copy()

    print("Adding unintegrated UMAP...")
    add_umap(input, counts=False)
    suffix_embeddings(input)

    output = run_symphony(input, seed)

    print("Adding integrated UMAP...")
    add_umap(output, use_rep="X_emb")

    print("Restoring counts matrix...")
    output.X = counts
    # Delete the log1p metadata so scanpy doesn't think we have log transformed data
    del output.uns["log1p"]

    print(f"Writing output to '{out_dir}'...")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    adata_file = join(out_dir, "adata.h5ad")
    print(f"Saving minimised AnnData to '{adata_file}'...")
    output_min = minimise_anndata(
        output, obs=["Batch", "Label", "Unseen"], obsm=["X_emb", "X_pca"], var=["mean", "std"], varm=["PCs"], uns=["Species", "harmony"]
    )
    output_min.write_h5ad(adata_file)

    umap_file = join(out_dir, "umap-unintegrated.png")
    umap = plot_embedding(output, basis="X_umap_unintegrated")
    print(f"Saving unintegrated UMAP plot to '{umap_file}'...")
    umap.savefig(umap_file, dpi=300, bbox_inches="tight")

    umap_file = join(out_dir, "umap-integrated.png")
    umap = plot_embedding(output)
    print(f"Saving integrated UMAP plot to '{umap_file}'...")
    umap.savefig(umap_file, dpi=300, bbox_inches="tight")

    print("Done!")


if __name__ == "__main__":
    main()
