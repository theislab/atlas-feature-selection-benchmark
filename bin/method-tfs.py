#!/usr/bin/env python

"""
Select transcription factor features

Usage:
    method-tfs.py --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --tfs-file=<path>    Path to TSV file containing transcription factors.
    --out-file=<path>    Path to output file.
"""


def select_tf_features(adata, tfs):
    """
    Select random features in a dataset

    Parameters
    ----------
    adata
        AnnData object
    tfs
        List of transcription factors

    Returns
    -------
    DataFrame containing the selected features
    """

    from pandas import DataFrame

    print("Selecting transcription factor features...")
    tfs = [gene for gene in tfs if gene in adata.var_names]
    selected_features = list(set(tfs))
    if len(selected_features) == 0:
        raise ValueError("No TF genes found in dataset")
    print(f"Selected {len(selected_features)} TF features")

    output = DataFrame(selected_features, columns=["Feature"])

    return output


def read_tf_genes(file, var_names, species):
    """
    Read transcription factor genes from a TSV file

    Parameters
    ----------
    file
        Path to TSV file containing transcription factor genes
    var_names
        List of the variable names for the dataset to select features for
    species
        The species of the data to be scored, either 'human' or 'mouse'

    Returns
    -------
    A list of transcription factor genes
    """

    from pandas import read_csv

    print(f"Reading TF genes from '{file}'...")
    tfs_df = read_csv(file, sep="\t")

    print(f"Selecting genes for species '{species}'...")
    if species.lower() == "human":
        tfs_df = tfs_df[tfs_df["Species"] == "Human"]
        ensembl_prefix = "ENSG"
    elif species.lower() == "mouse":
        tfs_df = tfs_df[tfs_df["Species"] == "Mouse"]
        ensembl_prefix = "ENSMUSG"
    else:
        from warnings import warn

        warn(
            f"'{species}' is not a valid species ('human', 'mouse'). Returning empty gene list."
        )
        return []

    is_ensembl = all(var_name.startswith(ensembl_prefix) for var_name in var_names)
    if is_ensembl:
        print(f"ENSEMBL var names detected, using ENSEMBL IDs")
        id_col = "ENSEMBL"
    else:
        print("ENSEMBL var names not found, using gene names")
        id_col = "Gene"

    tfs_df = tfs_df.dropna(subset=[id_col])
    tfs = list(tfs_df[id_col])

    return tfs

def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad

    args = docopt(__doc__)

    file = args["<file>"]
    tfs_file = args["--tfs-file"]
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    print("Read data:")
    print(input)
    tfs = read_tf_genes(tfs_file, input.var_names, input.uns["Species"])
    output = select_tf_features(input, tfs)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
