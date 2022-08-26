#!/usr/bin/env python

"""
Map a query dataset to an scVI reference model

Usage:
    map-scvi.py --reference=<path> --out-dir=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --reference=<path>   Path to reference model directory.
    --out-dir=<path>     Path to output directory.
"""

import scvi


def map_query_scVI(reference, query):
    """
    Map a query to an scVI reference model

    Parameters
    ----------
    reference
        The reference model to map to
    query
        AnnData object containing the query dataset

    Returns
    -------
    Mapped query model
    """

    print("Creating scVI query model...")
    model = scvi.model.SCVI.load_query_data(
        query,
        reference,
    )
    print(model)
    model.view_anndata_setup()

    print("Training scVI query model...")
    model.train(max_epochs=200, plan_kwargs=dict(weight_decay=0.0))
    print(model)

    return model


def add_integrated_embeddings(model, adata):
    """
    Add embeddings from an integration model to an AnnData

    Parameters
    ----------
    model
        The trained integration model
    adata
        The AnnData to add embeddings to

    Returns
    -------
    AnnData containing integrated embeddings
    """

    print("Storing unintegrated embeddings...")
    for key in adata.obsm_keys():
        print(f"Storing {key}...")
        adata.obsm[key + "_unintegrated"] = adata.obsm[key].copy()

    print("Adding integrated embedding...")
    adata.obsm["X_scVI"] = model.get_latent_representation(adata)

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
        basis=basis,
        color=["Dataset", "Batch", "Label"],
        legend_fontsize="small",
        legend_fontweight="light",
        title=["Dataset", "Batch", "Label"],
        add_outline=True,
        outline_width=(0.1, 0.05),
        ncols=1,
        show=False,
        return_fig=True,
    )

    return plt


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad
    from os.path import join

    args = docopt(__doc__)

    file = args["<file>"]
    reference_dir = args["--reference"]
    out_dir = args["--out-dir"]

    print(f"Reading reference model from '{reference_dir}'...")
    reference = scvi.model.SCVI.load(reference_dir)
    del reference.adata.obsm  # Clear the reference embeddings
    print("Read model:")
    print(reference)
    print(f"Reading query from '{file}'...")
    input = read_h5ad(file)
    print("Read query:")
    print(input)
    print("Subsetting query to selected features...")
    input = input[:, reference.adata.var_names].copy()
    print(input)
    output = map_query_scVI(reference, input)
    print("Concatenating reference and query...")
    input.obs["Dataset"] = "Query"
    reference.adata.obs["Dataset"] = "Reference"
    full = input.concatenate(reference.adata)
    print("Adding unintegrated UMAP...")
    add_umap(full)
    add_integrated_embeddings(output, full)
    print(f"Writing output to '{out_dir}'...")
    output.save(out_dir, save_anndata=True, overwrite=True)
    umap_file = join(out_dir, "umap-unintegrated.png")
    umap = plot_umap(full, basis="X_umap_unintegrated")
    print(f"Saving unintegrated UMAP plot to '{umap_file}'...")
    umap.savefig(umap_file, dpi=300, bbox_inches="tight")
    umap_file = join(out_dir, "umap-integrated.png")
    umap = plot_umap(full)
    print(f"Saving integrated UMAP plot to '{umap_file}'...")
    umap.savefig(umap_file, dpi=300, bbox_inches="tight")
    print("Done!")


if __name__ == "__main__":
    main()
