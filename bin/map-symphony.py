#!/usr/bin/env python

"""
Map a query dataset to a Symphony reference model

Usage:
    map-symphony.py--reference=<file> --reference-model=<path> --out-dir=<path> [options] <file>

Options:
    -h --help                  Show this screen.
    --reference=<file>         Path to reference AnnData file.
    --reference-model=<path>   Path to reference model directory.
    --out-dir=<path>           Path to output directory.
"""

from functions.integration import (
    add_umap,
    suffix_embeddings,
    plot_embedding,
)


def map_query_symphony(reference, query):
    """
    Map a query to a Symphony reference model

    Parameters
    ----------
    reference
        AnnData object containing the reference to map to
    query
        AnnData object containing the query dataset

    Returns
    -------
    AnnData with mapped query dataset
    """

    from symphonypy.tools import map_embedding

    print("Mapping query using Symphony...")
    map_embedding(
        query,
        reference,
        key="Batch",
        use_genes_column=None,
        transferred_adjusted_basis="X_emb"
    )

    return query
def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad
    from scanpy.preprocessing import normalize_total, log1p
    from os.path import join
    from functions.anndata import minimise_anndata
    import os

    args = docopt(__doc__)

    file = args["<file>"]
    reference_file = args["--reference"]
    reference_dir = args["--reference-model"]
    out_dir = args["--out-dir"]

    print(f"Reading reference AnnData from '{reference_file}'...")
    reference_adata = read_h5ad(reference_file)

    print(f"Reading reference model from '{reference_dir}'...")
    model_adata = read_h5ad(join(reference_dir, "adata.h5ad"))
    print("Normalising reference...")
    normalize_total(reference_adata, target_sum=1e4)
    log1p(reference_adata)
    model_adata.X = reference_adata[:, model_adata.var_names].X.copy()
    del reference_adata
    print("Read model:")
    print(model_adata)

    print(f"Reading query from '{file}'...")
    input = read_h5ad(file)
    print("Storing query counts matrix...")
    counts = input.X.copy()
    print("Normalising query...")
    normalize_total(input, target_sum=1e4)
    log1p(input)
    print("Read query:")
    print(input)

    print("Subsetting query to selected features...")
    query = input[:, model_adata.var_names].copy()
    counts = counts[:, input.var_names.isin(model_adata.var_names)].copy()
    del input
    print(query)

    print("Adding unintegrated UMAP to query..")
    add_umap(query, counts=False)
    suffix_embeddings(query)

    output = map_query_symphony(model_adata, query)

    print("Adding integrated UMAP to query...")
    add_umap(output, use_rep="X_emb")

    print("Concatenating reference and query...")
    output.obs["Dataset"] = "Query"
    model_adata.obs["Dataset"] = "Reference"
    full = output.concatenate(model_adata)

    print("Adding unintegrated UMAP...")
    add_umap(full, counts=False)
    suffix_embeddings(full)
    print("Adding integrated UMAP...")
    full.obsm["X_emb"] = full.obsm["X_emb_unintegrated"] # Restore the name we just replaced
    add_umap(full, use_rep="X_emb")

    print("Restoring query counts matrix...")
    output.X = counts
    # Delete the log1p metadata so scanpy doesn't think we have log transformed data
    del output.uns["log1p"]

    print(f"Writing output to '{out_dir}'...")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    adata_file = join(out_dir, "adata.h5ad")
    print(f"Saving minimised AnnData to '{adata_file}'...")
    output_min = minimise_anndata(
        output, obs=["Batch", "Label", "Unseen"], obsm=["X_emb"], uns=["Species"]
    )
    output_min.write_h5ad(adata_file)

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
