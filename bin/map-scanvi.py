#!/usr/bin/env python

"""
Map a query dataset to an scANVI reference model

Usage:
    map-scanvi.py --reference=<path> --out-dir=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --reference=<path>   Path to reference model directory.
    --out-dir=<path>     Path to output directory.
"""

import scvi
from functions.integration import (
    add_umap,
    add_integrated_embeddings,
    suffix_embeddings,
    plot_embedding,
)


def map_query_scANVI(reference, query):
    """
    Map a query to an scANVI reference model

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

    print("Creating scANVI query model...")
    model = scvi.model.SCANVI.load_query_data(
        query,
        reference,
    )
    print(model)
    model.view_anndata_setup()

    print("Training scVI query model...")
    model.train(
        max_epochs=100,
        plan_kwargs=dict(weight_decay=0.0),
        check_val_every_n_epoch=10,
    )
    print(model)

    return model


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad
    from os.path import join

    args = docopt(__doc__)

    file = args["<file>"]
    reference_dir = args["--reference"]
    adata_file = join(reference_dir, "adata.h5ad")
    out_dir = args["--out-dir"]

    print(f"Reading reference AnnData from '{adata_file}'...")
    reference_adata = read_h5ad(adata_file)
    del reference_adata.obsm  # Clear the reference embeddings
    print("Subsetting reference to selected features...")
    model_adata = reference_adata[:, reference_adata.var["Selected"]].copy()

    print(f"Reading reference model from '{reference_dir}'...")
    reference = scvi.model.SCANVI.load(reference_dir, adata=model_adata)
    print("Read model:")
    print(reference)
    print(f"Reading query from '{file}'...")
    input = read_h5ad(file)
    print("Read query:")
    print(input)
    print("Subsetting query to selected features...")
    query = input[:, model_adata.var_names].copy()
    print(query)

    output = map_query_scANVI(reference, query)

    print("Adding embeddings to query...")
    add_umap(output.adata)
    suffix_embeddings(output.adata)
    add_integrated_embeddings(output, output.adata)

    print("Concatenating reference and query...")
    output.adata.obs["Dataset"] = "Query"
    model_adata.obs["Dataset"] = "Reference"
    full = output.adata.concatenate(model_adata)

    print("Adding unintegrated UMAP...")
    add_umap(full)
    suffix_embeddings(full)
    add_integrated_embeddings(output, full)

    print(f"Writing output to '{out_dir}'...")
    input.obsm = full[full.obs["Dataset"] == "Query"].obsm.copy()
    output.save(out_dir, save_anndata=False, overwrite=True)
    input.write_h5ad(join(out_dir, "adata.h5ad"))
    # Set unseen labels to string for plotting
    full.obs["Unseen"] = full.obs["Unseen"].astype(str)
    umap_file = join(out_dir, "umap-unintegrated.png")
    umap = plot_embedding(
        full,
        basis="X_umap_unintegrated",
        groups=["Dataset", "Batch", "Label", "Unseen"],
    )
    print(f"Saving unintegrated UMAP plot to '{umap_file}'...")
    umap.savefig(umap_file, dpi=300, bbox_inches="tight")
    umap_file = join(out_dir, "umap-integrated.png")
    umap = plot_embedding(full, groups=["Dataset", "Batch", "Label", "Unseen"])
    print(f"Saving integrated UMAP plot to '{umap_file}'...")
    umap.savefig(umap_file, dpi=300, bbox_inches="tight")
    print("Done!")


if __name__ == "__main__":
    main()
