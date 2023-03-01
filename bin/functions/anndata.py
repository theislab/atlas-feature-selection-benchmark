def minimise_anndata(
    adata, X=False, obs=None, layers=None, var=None, uns=None, obsm=None, varm=None
):
    """
    Minimise an AnnData object by removing unnecessary data

    Parameters
    ----------
    adata
        AnnData object to minimise
    X
        Whether to keep the X matrix
    layers
        List of keys in adata.layers to keep
    obs
        List of columns in adata.obs to keep
    var
        List of columns in adata.var to keep
    uns
        List of keys in adata.uns to keep
    obsm
        List of keys in adata.obsm to keep
    varm
        List of keys in adata.varm to keep

    Returns
    -------
    Minimised AnnData object
    """

    from anndata import AnnData

    adata_min = AnnData(shape=adata.shape)
    adata_min.var_names = adata.var_names.copy()
    adata_min.obs_names = adata.obs_names.copy()

    if X:
        adata_min.X = adata.X.copy()

    if layers is not None:
        for key in layers:
            adata_min.layers[key] = adata.layers[key].copy()

    if obs is not None:
        for key in obs:
            adata_min.obs[key] = adata.obs[key]

    if var is not None:
        for key in var:
            adata_min.var[key] = adata.var[key]

    if uns is not None:
        for key in uns:
            adata_min.uns[key] = adata.uns[key]

    if obsm is not None:
        for key in obsm:
            adata_min.obsm[key] = adata.obsm[key]

    if varm is not None:
        for key in varm:
            adata_min.varm[key] = adata.varm[key]

    adata_min = adata_min.copy()

    return adata_min
