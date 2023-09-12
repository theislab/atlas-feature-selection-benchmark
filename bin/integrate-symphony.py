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

    from scanpy.preprocessing import normalize_total, log1p, scale
    from scanpy.tools import pca
    from symphonypy.preprocessing import harmony_integrate

    print("Storing counts matrix...")
    counts = adata.X.copy()

    print("Normalising...")
    normalize_total(adata, target_sum=1e4)

    print("Log transforming...")
    log1p(adata)

    print("Scaling...")
    scale(adata, max_value=10)

    print("Calculating PCA...")
    pca(adata, n_comps=30, zero_center=False)

    print("Integrating with Harmony...")
    print(f"Using random seed: {seed}")
    harmony_integrate(
        adata, key="Batch", ref_basis_adjusted="X_emb", verbose=True, random_seed=seed
    )

    print("Restoring counts matrix...")
    adata.X = counts

    return adata


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad
    from pandas import read_csv
    from os.path import join
    from functions.anndata import minimise_anndata
    from pickle import dump
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

    print(f"Reading selected features from '{features_file}'...")
    features = read_csv(features_file, sep="\t")
    print("Read features:")
    print(features)

    print(f"Subsetting to {features.shape[0]} selected features...")
    input = input[:, input.var_names.isin(features["Feature"])].copy()

    print("Adding unintegrated UMAP...")
    add_umap(input)
    suffix_embeddings(input)

    output = run_symphony(input, seed)

    print("Adding integrated UMAP...")
    add_umap(output, use_rep="X_emb")

    print(f"Writing output to '{out_dir}'...")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    harmony_model = output.uns["harmony"].copy()
    with open(join(out_dir, "harmony.pkl"), "wb") as harmony_file:
        print(f"Saving Harmony model to '{harmony_file}'...")
        dump(harmony_model, harmony_file)

    adata_file = join(out_dir, "adata.h5ad")
    print(f"Saving minimised AnnData to '{adata_file}'...")
    output_min = minimise_anndata(
        output, obs=["Batch", "Label", "Unseen"], obsm=["X_emb"], uns=["Species"]
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
