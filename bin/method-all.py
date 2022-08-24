#!/usr/bin/env python

"""
Select all features

Usage:
    method-all.py --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
"""


def select_all_features(adata):
    """
    Select all the features in a dataset

    Parameters
    ----------
    adata
        AnnData object

    Returns
    -------
    DataFrame containing the selected features
    """

    from pandas import DataFrame

    print("Selecting all features...")
    selected_features = DataFrame({"Feature" : adata.var_names})

    return selected_features


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad

    args = docopt(__doc__)

    file = args["<file>"]
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    print("Read data:")
    print(input)
    output = select_all_features(input)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
