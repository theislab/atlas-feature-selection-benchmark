#!/usr/bin/env python

"""
Integrate a dataset using scANVI

Usage:
    integrate-scanvi.py --out-dir=<path> [options] <dir>

Options:
    -h --help            Show this screen.
    --out-dir=<path>     Path to output directory.
"""

import scvi
from functions import add_integrated_embeddings, plot_embedding


def run_scANVI(scvi_model):
    """
    Integrate a dataset using scANVI

    Parameters
    ----------
    scvi_model
        A trained scVI model

    Returns
    -------
    Integrated scANVI model
    """

    print("Setting reference labels...")
    adata = scvi_model.adata
    adata.obs["ReferenceLabel"] = adata.obs["Label"].values

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

    args = docopt(__doc__)

    dir = args["<dir>"]
    out_dir = args["--out-dir"]

    print(f"Reading model from '{dir}'...")
    input = scvi.model.SCVI.load(dir)
    print("Read model:")
    print(input)
    output = run_scANVI(input)
    print("Storing scVI embeddings...")
    output.adata.obsm["X_scVI"] = output.adata.obsm["X_emb"].copy()
    output.adata.obsm["X_umap_scVI"] = output.adata.obsm["X_umap"].copy()
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
