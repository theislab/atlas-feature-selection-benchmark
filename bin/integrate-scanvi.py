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

def run_scANVI(scvi_model):
    """
    A function that performs analysis

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

def add_integrated_embeddings(model):
    """
    Add embeddings from an integration model to an AnnData

    Parameters
    ----------
    model
        The trained integration model containing an AnnData object

    Returns
    -------
    Model with AnnData containing integrated embeddings
    """

    adata = model.adata

    print("Storing unintegrated embeddings...")
    for key in adata.obsm_keys():
        print(f"Storing {key}...")
        adata.obsm[key + "_unintegrated"] = adata.obsm[key].copy()

    print("Adding integrated embedding...")
    adata.obsm["X_scVI"] = model.get_latent_representation()

    add_umap(adata, use_rep="X_scVI")

    return None

def add_umap(adata, use_rep=None):
    """
    Add a UMAP embedding to an AnnData

    Parameters
    ----------
    adata
        AnnData object to add UMAP to
    use_rep
        The base representation to use. If `None` then PCA is calculated first.

    Returns
    -------
    Model with AnnData containing integrated embeddings
    """

    from scanpy.preprocessing import neighbors
    from scanpy.tools import pca, umap

    if use_rep is None:
        print("Calculating PCA...")
        pca(adata)
        use_rep = "X_pca"

    print(f"Calculating nearest neighbours using '{use_rep}'...")
    neighbors(adata, use_rep=use_rep)

    print("Calculating UMAP...")
    umap(adata)

    return None

def plot_umap(adata, basis="X_umap"):
    """
    Plot a UMAP

    Parameters
    ----------
    adata
        AnnData object containing the UMAP to plot

    Returns
    -------
    matplotlib object
    """

    from scanpy.plotting import embedding

    print(f"Plotting UMAP using '{basis}'...")
    plt = embedding(
        adata,
        basis = basis,
        color=["Batch", "Label"],
        legend_fontsize="small",
        legend_fontweight="light",
        title=["Batch", "Label"],
        add_outline=True,
        outline_width=(0.1, 0.05),
        ncols=1,
        show=False,
        return_fig=True
    )

    return plt

def main():
    """The main script function"""
    from docopt import docopt
    from os.path import join

    args = docopt(__doc__)

    dir = args["<dir>"]
    out_dir = args["--out-dir"]

    print(f"Reading model from '{dir}'...")
    input = scvi.model.SCVI.load(dir)
    del input.adata.obsm # Clear the embeddings from scVI
    print("Read model:")
    print(input)
    output = run_scANVI(input)
    print("Adding unintegrated UMAP...")
    add_umap(output.adata)
    add_integrated_embeddings(output)
    print(f"Writing output to '{out_dir}'...")
    output.save(out_dir, save_anndata=True, overwrite=True)
    umap_file = join(out_dir, "umap-unintegrated.png")
    umap = plot_umap(output.adata, basis="X_umap_unintegrated")
    print(f"Saving unintegrated UMAP plot to '{umap_file}'...")
    umap.savefig(umap_file, dpi=300, bbox_inches="tight")
    umap_file = join(out_dir, "umap-integrated.png")
    umap = plot_umap(output.adata)
    print(f"Saving integrated UMAP plot to '{umap_file}'...")
    umap.savefig(umap_file, dpi=300, bbox_inches="tight")
    print("Done!")


if __name__ == "__main__":
    main()
