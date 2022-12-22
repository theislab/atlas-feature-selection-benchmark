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
    --species=<str>            Species of the dataset.
    --reference-out=<path>     Path to reference output file.
    --query-out=<path>         Path to query output file.
"""


def prepare_dataset(adata_raw, name, batch_col, label_col, query_batches, species):
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
    species
        Species of the dataset

    Returns
    -------
    A tuple containing the prepared reference and query AnnData objects
    """

    from scanpy.preprocessing import filter_cells, filter_genes
    from anndata import AnnData

    print(f"Preparing '{name}' dataset...")
    print_adata(adata_raw, "RAW DATA", batch_col, label_col)

    print("Creating new AnnData...")
    adata = AnnData(X=adata_raw.X.copy())
    adata.obs_names = adata_raw.obs_names.copy()
    adata.var_names = adata_raw.var_names.copy()
    adata.obs["Batch"] = adata_raw.obs[batch_col]
    adata.obs["Label"] = adata_raw.obs[label_col]
    adata.uns["Species"] = species
    del adata_raw

    print("Removing cells with less than 100 counts...")
    filter_cells(adata, min_counts=100)

    print("Removing cells with less than 100 genes...")
    filter_cells(adata, min_genes=100)

    print("Splitting reference and query...")
    print(
        f"Selecting {len(query_batches)} batch(es) as the query: {', '.join(query_batches)}"
    )
    is_query = adata.obs["Batch"].isin(query_batches)
    reference = adata[~is_query].copy()
    reference.obs["Batch"] = reference.obs["Batch"].cat.remove_unused_categories()
    query = adata[is_query].copy()
    query.obs["Batch"] = query.obs["Batch"].cat.remove_unused_categories()
    del adata

    print("Removing labels with fewer than 20 cells...")
    label_counts = reference.obs["Label"].value_counts()
    keep_labels = list(label_counts.index[label_counts >= 20])
    reference = reference[reference.obs["Label"].isin(keep_labels)].copy()
    reference.obs["Label"] = reference.obs["Label"].cat.remove_unused_categories()
    label_counts = query.obs["Label"].value_counts()
    # Remove uncommon labels in the query but
    # make sure we keep can labels from the reference
    keep_labels = set(keep_labels + list(label_counts.index[label_counts >= 20]))
    query = query[query.obs["Label"].isin(keep_labels)].copy()
    query.obs["Label"] = query.obs["Label"].cat.remove_unused_categories()

    print("Removing genes with 0 counts in the reference...")
    filter_genes(reference, min_counts=1)
    query = query[:, reference.var_names].copy()

    print("Removing QC stats...")
    reference.obs = reference.obs.drop(["n_counts", "n_genes"], axis=1)
    reference.var = reference.var.drop(["n_counts"], axis=1)
    query.obs = query.obs.drop(["n_counts", "n_genes"], axis=1)

    print_adata(reference, "REFERENCE DATA", "Batch", "Label")
    print_adata(query, "QUERY DATA", "Batch", "Label")

    return (reference, query)


def print_adata(adata, title, batch_col, label_col):
    """
    Print an AnnData object with information about batches, labels and species

    Parameters
    ----------
    adata
        AnnData object to print
    title
        Title for the object
    batch_col
        Name of the column of adata.obs containing batch information
    label_col
        Name of the column of adata.obs containing label information
    """

    print()
    print(f"====== {title} =======")
    print(f"Cells: {adata.n_obs}")
    print(f"Genes: {adata.n_vars}")
    print("-----------------------")
    print(f"Batches: '{batch_col}' ({len(adata.obs[batch_col].cat.categories)})")
    print(adata.obs[batch_col].value_counts().to_string())
    print("-----------------------")
    print(f"Labels: '{label_col}' ({len(adata.obs[label_col].cat.categories)})")
    print(adata.obs[label_col].value_counts().to_string())
    print("-----------------------")
    if "Species" in adata.uns_keys():
        print(f"Species: '{adata.uns['Species']}'")
    else:
        print("Species not set")
    print("-----------------------")
    print("Object:")
    print(adata)
    print("=======================")
    print()


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
    species = args["--species"]
    reference_out = args["--reference-out"]
    query_out = args["--query-out"]

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    reference, query = prepare_dataset(
        input, name, batch_col, label_col, query_batches, species
    )
    print(f"Writing reference to '{reference_out}'...")
    reference.write_h5ad(reference_out)
    print(f"Writing query to '{query_out}'...")
    query.write_h5ad(query_out)
    print("Done!")


if __name__ == "__main__":
    main()
