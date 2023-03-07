#!/usr/bin/env python

"""
Integrate a dataset using scANVI

Usage:
    integrate-scanvi.py --scvi=<dir> --out-dir=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --scvi=<path>        Path to scVI model output directory.
    --out-dir=<path>     Path to output directory.
    --seed=<int>         Random seed to use [default: 0].
"""

import scvi
from functions.integration import (
    add_umap,
    add_integrated_embeddings,
    suffix_embeddings,
    plot_embedding,
)


def run_scANVI(scvi_model, seed):
    """
    Integrate a dataset using scANVI

    Parameters
    ----------
    scvi_model
        A trained scVI model
    seed
        Random seed to use

    Returns
    -------
    Integrated scANVI model
    """

    print(f"Setting random seed to {seed}...")
    scvi.settings.seed = seed

    print("Creating scANVI model...")
    model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        unlabeled_category="Unknown",
        labels_key="ReferenceLabel",
    )
    print(model)
    model.view_anndata_setup()

    print("Training scANVI model...")
    model.train(max_epochs=20, n_samples_per_label=100)
    print(model)

    return model


def main():
    """The main script function"""
    from docopt import docopt
    from os.path import join
    from functions.anndata import minimise_anndata

    args = docopt(__doc__)

    file = args["<file>"]
    scvi_dir = args["--scvi"]
    out_dir = args["--out-dir"]
    seed = int(args["--seed"])

    print(f"Reading AnnData from '{file}'...")
    adata = scvi.data.read_h5ad(file)
    print(adata)

    print(f"Reading model from '{scvi_dir}'...")
    model_adata = scvi.data.read_h5ad(join(scvi_dir, "adata.h5ad"))
    model_adata.X = adata[:, model_adata.var_names].X.copy()
    model_adata.obs["ReferenceLabel"] = model_adata.obs[
        "Label"
    ].cat.remove_unused_categories()
    input = scvi.model.SCVI.load(scvi_dir, adata=model_adata)
    print("Read model:")
    print(input)

    output = run_scANVI(input, seed)

    print("Storing scVI embedding...")
    scvi_emb = output.adata.obsm["X_emb"].copy()
    del output.adata.obsm["X_emb"]

    print("Adding unintegrated UMAP...")
    add_umap(output.adata)
    suffix_embeddings(output.adata)
    add_integrated_embeddings(output, output.adata)

    print(f"Writing output to '{out_dir}'...")
    output.save(out_dir, save_anndata=False, overwrite=True)
    output.adata.obsm["X_scVI"] = scvi_emb
    output_min = minimise_anndata(
        output.adata,
        obs=["Batch", "Label", "Unseen", "ReferenceLabel"],
        obsm=["X_emb", "X_scVI"],
        uns=["Species"],
    )
    output_min.write_h5ad(join(out_dir, "adata.h5ad"))

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
