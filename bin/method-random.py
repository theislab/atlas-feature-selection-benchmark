#!/usr/bin/env python

"""
Select random features

Usage:
    method-random.py --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
    --n-features=<int>   Number of features to select [default: 1000].
    --seed=<int>         Random seed to use [default: 1].
"""


def select_random_features(adata, n, seed):
    """
    Select random features in a dataset

    Parameters
    ----------
    adata
        AnnData object
    n
        Number of features to select

    Returns
    -------
    DataFrame containing the selected features
    """

    from pandas import DataFrame

    if adata.n_vars < n:
        import warnings

        warnings.warn(
            "Number of features to select is greater than the number present, setting n to adata.n_vars"
        )
        n = adata.n_vars

    print(f"Selecting {n} random features with seed {seed}...")
    selected_features = DataFrame({"Feature": adata.var_names})
    selected_features = selected_features.sample(n=n, random_state=seed)

    return selected_features


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad

    args = docopt(__doc__)

    file = args["<file>"]
    n_features = int(args["--n-features"])
    seed = int(args["--seed"])
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    print("Read data:")
    print(input)
    output = select_random_features(input, n_features, seed)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
