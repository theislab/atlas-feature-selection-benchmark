#!/usr/bin/env python

"""
Prepare a dataset for the benchmarking pipeline.

Usage:
    prepare-dataset.py --name=<str> --batch-col=<str> --label-col=<str> --query-batches=<str> --reference-out=<path> --query-out=path [options] <file>

Options:
    -h --help                  Show this screen.
    --name=<str>               Name of the dataset.
    --batch-col=<str>          Name of the obs column containing batch information.
    --label-col=<str>          Name of the obs column containing label information.
    --query-batches=<str>      Comma separated list of batches to be used as the query.
    --reference-out=<path>     Path to reference output file.
    --query-out=<path>         Path to query output file.
"""


def prepare_dataset(adata_raw, name, batch_col, label_col, query_batches):
    """
    Prepare a dataset for the benchmarking pipeline

    Parameters
    ----------
    adata_raw
        AnnData object containing the raw dataset to prepare
    name
        Name of the dataset
    batch_col
        Name of the column of adata.obs containing batch information
    label_col
        Name of the column of adata.obs containing label information
    query_batches
        List containing the batch names for the query

    Returns
    -------
    A tuple containing the prepared reference and query AnnData objects
    """

    from scanpy.preprocessing import filter_cells, filter_genes
    from anndata import AnnData

    print(f"Preparing {name} dataset...")
    print()
    print("====== RAW DATA =======")
    print(f"Cells: {adata_raw.n_obs}")
    print(f"Genes: {adata_raw.n_vars}")
    print(
        f"Batches ({len(adata_raw.obs[batch_col].cat.categories)}): {batch_col} ({', '.join(adata_raw.obs[batch_col].cat.categories)})"
    )
    print(
        f"Labels ({len(adata_raw.obs[label_col].cat.categories)}): {label_col} ({', '.join(adata_raw.obs[label_col].cat.categories)})"
    )
    print("Object:")
    print(adata_raw)
    print("=======================")
    print()

    print("Creating new AnnData...")
    adata = AnnData(X=adata_raw.X.copy())
    adata.obs_names = adata_raw.obs_names.copy()
    adata.var_names = adata_raw.var_names.copy()
    adata.obs["Batch"] = adata_raw.obs[batch_col]
    adata.obs["Label"] = adata_raw.obs[label_col]
    del adata_raw

    print("Removing cells with less than 100 counts...")
    filter_cells(adata, min_counts=100)

    print("Removing cells with less than 100 genes...")
    filter_cells(adata, min_genes=100)

    print("Removing genes with 0 counts...")
    filter_genes(adata, min_counts=1)

    print("Removing unused labels...")
    adata.obs["Label"] = adata.obs["Label"].cat.remove_unused_categories()

    print("Removing QC stats...")
    adata.obs = adata.obs.drop(["n_counts", "n_genes"], axis=1)
    adata.var = adata.var.drop(["n_counts"], axis=1)

    print("Splitting reference and query...")
    print(
        f"Selecting {len(query_batches)} batches as the query: {', '.join(query_batches)}"
    )
    is_query = adata.obs["Batch"].isin(query_batches)
    adata[adata.obs["Batch"].isin(query_batches)]
    reference = adata[~is_query].copy()
    reference.obs["Batch"] = reference.obs["Batch"].cat.remove_unused_categories()
    query = adata[is_query].copy()
    query.obs["Batch"] = query.obs["Batch"].cat.remove_unused_categories()

    print()
    print("====== REFERENCE DATA =======")
    print(f"Cells: {reference.n_obs}")
    print(f"Genes: {reference.n_vars}")
    print(
        f"Batches ({len(reference.obs['Batch'].cat.categories)}): Batch ({', '.join(reference.obs['Batch'].cat.categories)})"
    )
    print(
        f"Labels ({len(reference.obs['Label'].cat.categories)}): Label ({', '.join(reference.obs['Label'].cat.categories)})"
    )
    print("Object:")
    print(reference)
    print("============================")
    print()
    print("====== QUERY DATA =======")
    print(f"Cells: {query.n_obs}")
    print(f"Genes: {query.n_vars}")
    print(
        f"Batches ({len(query.obs['Batch'].cat.categories)}): Batch ({', '.join(query.obs['Batch'].cat.categories)})"
    )
    print(
        f"Labels ({len(query.obs['Label'].cat.categories)}): Label ({', '.join(query.obs['Label'].cat.categories)})"
    )
    print("Object:")
    print(query)
    print("=========================")
    print()

    return (reference, query)


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad

    args = docopt(__doc__)

    file = args["<file>"]
    name = args["--name"]
    batch_col = args["--batch-col"]
    label_col = args["--label-col"]
    query_batches = args["--query-batches"].split(",")
    reference_out = args["--reference-out"]
    query_out = args["--query-out"]

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    reference, query = prepare_dataset(input, name, batch_col, label_col, query_batches)
    print(f"Writing reference to '{reference_out}'...")
    reference.write_h5ad(reference_out)
    print(f"Writing query to '{query_out}'...")
    query.write_h5ad(query_out)
    print("Done!")


if __name__ == "__main__":
    main()
