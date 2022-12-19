#!/usr/bin/env python

"""
Select features using the anti-correlation method

Usage:
    method-anticor.py --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
"""


def select_anticor_features(input, temp_dir):
    """
    Select features base on anti-correlations

    Parameters
    ----------
    input
        AnnData to select features from
    temp_dir
        Path to a temporary directory to use for calculations

    Returns
    -------
    DataFrame containing the selected features
    """

    from anticor_features.anticor_features import get_anti_cor_genes
    from os import makedirs
    from scanpy.preprocessing import log1p, normalize_total

    print("Normalising expression...")
    normalize_total(input)
    log1p(input)

    print("Selecting anticor features...")
    print(f"Using temporary directory: '{temp_dir}'")
    makedirs(temp_dir, exist_ok=True)
    anticor = get_anti_cor_genes(
        input.X.T,
        input.var.index.tolist(),
        pre_remove_pathways=[],    # Turn off pathway filtering
        scratch_dir=temp_dir       # Set temporary directory
    )

    print("Selecting top features...")
    selected = anticor[anticor["selected"]==True].sort_values(by="FDR").copy()
    selected["Feature"] = selected["gene"]

    return selected


def main():
    """The main script function"""
    from docopt import docopt
    from anndata import read_h5ad
    from tempfile import gettempdir
    from os.path import basename, join, splitext

    args = docopt(__doc__)

    file = args["<file>"]
    out_file = args["--out-file"]

    temp_dir = join(gettempdir(), "anticor", splitext(basename(file))[0])

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    print("Read data:")
    print(input)
    output = select_anticor_features(input, temp_dir)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
