def add_umap(adata, use_rep=None, counts=True):
    """
    Add a UMAP embedding to an AnnData

    Parameters
    ----------
    adata
        AnnData object to add UMAP to
    use_rep
        The base representation to use. If `None` then PCA is calculated first.
    counts
        If `use_rep=None` then whether `adata.X` contains raw counts or not.

    Returns
    -------
    Model with AnnData containing integrated embeddings
    """

    from scanpy.preprocessing import neighbors
    from scanpy.tools import pca, umap

    if use_rep is None:
        if counts:
            print("Calculating normalised matrix for PCA...")
            from scanpy.preprocessing import normalize_total, log1p

            counts_mat = adata.X.copy()
            normalize_total(adata, target_sum=1e4)
            log1p(adata)

            print("Calculating PCA...")
            pca(adata)

            print("Restoring counts matrix...")
            adata.X = counts_mat
            # Delete the log1p metadata so scanpy doesn't think we have log transformed data
            del adata.uns["log1p"]

        else:
            print("Calculating PCA...")
            pca(adata)

        use_rep = "X_pca"

    print(f"Calculating nearest neighbours using '{use_rep}'...")
    neighbors(adata, use_rep=use_rep)

    print("Calculating UMAP...")
    umap(adata)

    return None


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
    AnnData containing integrated embeddings in adata.obsm["X_emb"] and
    adata.obsm["X_umap"]
    """

    print("Adding integrated embedding...")
    adata.obsm["X_emb"] = model.get_latent_representation(adata)

    add_umap(adata, use_rep="X_emb")

    return None


def suffix_embeddings(adata, suffix="_unintegrated"):
    print(f"Adding '{suffix}' suffix to embeddings...")
    for key in adata.obsm_keys():
        print(f"Storing '{key}' as '{key + suffix}'...")
        adata.obsm[key + suffix] = adata.obsm[key].copy()
        del adata.obsm[key]

    return None


def plot_embedding(adata, basis="X_umap", groups=["Batch", "Label"]):
    """
    Plot a UMAP

    Parameters
    ----------
    adata
        AnnData object containing the UMAP to plot
    basis
        The name of the embedding in adata.obsm to plot
    groups
        List of the columns in adata.obs to colour points by.

    Returns
    -------
    matplotlib object
    """

    from scanpy.plotting import embedding

    print(f"Plotting UMAP using '{basis}'...")
    plt = embedding(
        adata,
        basis=basis,
        color=groups,
        legend_fontsize="small",
        legend_fontweight="light",
        title=groups,
        add_outline=True,
        outline_width=(0.1, 0.05),
        ncols=1,
        show=False,
        return_fig=True,
    )

    return plt
