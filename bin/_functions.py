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

    print(f"Adding {suffix} suffix to embeddings...")
    for key in adata.obsm_keys():
        print(f"Storing {key}...")
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


def format_metric_results(dataset, method, integration, metric_type, metric, value):
    """
    Format metric results

    Parameters
    ----------
    dataset
        The name of the dataset the metric has been calculated for
    method
        The name of the method the metric has been calculated for
    integration
        The name of the integration the metric has been calculated for
    metric_type
        The type of the metric, either "Integration" or "Classification"
    metric
        The name of the metric that has been calculated
    value
        The value of the calculated metric

    Returns
    -------
    DataFrame containing the formatted results
    """

    from pandas import DataFrame

    if not metric_type in ["Integration", "Classification", "Mapping"]:
        raise ValueError(
            "'metric_type' must be one of 'Integration', 'Classification', or 'Mapping'"
        )

    if value < 0 or value > 1:
        raise ValueError("'score' must be between 0 and 1")

    print("Formatting metric results...")
    results = DataFrame(
        [
            {
                "Dataset": dataset,
                "Method": method,
                "Integration": integration,
                "Type": metric_type,
                "Metric": metric,
                "Value": value,
            }
        ]
    )

    return results
