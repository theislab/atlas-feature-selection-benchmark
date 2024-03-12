#!/usr/bin/env python

"""
Combine a set of feature selection output files. If only one file is given it is
assumed that is a directory to search for TSV feature selection output files.

Usage:
    combine-features.py --out-file=<path> [options] <file>...

Options:
    -h --help                Show this screen.
    --out-file=<path>        Path to combined metrics output file.
"""


def combine_features(files):
    """
    Combine a set of features output files

    Parameters
    ----------
    files
        List of paths to feature output files

    Returns
    -------
    DataFrame indicating which features are present in which files
    """

    from pandas import read_csv, DataFrame
    import os

    total_files = len(files)
    progress_step = total_files // 10

    # Handle cases with less than 10 files
    if progress_step == 0:
        progress_step = 1

    print(f"Reading features from {total_files} files...")

    file_names = [os.path.basename(path) for path in files]

    # Dictionary to store features and their presence across files
    feature_presence = {}

    for i, file_path in enumerate(files, 1):

        file_df = read_csv(file_path, sep="\t")
        file_features = set(file_df["Feature"])
        file_name = os.path.basename(file_path)

        for feature in file_features:
            if feature not in feature_presence:
                feature_presence[feature] = [False] * total_files
            feature_presence[feature][file_names.index(file_name)] = True

        if i % progress_step == 0 or i == total_files:
            print(
                f"Progress: {i}/{total_files} files processed ({i/total_files * 100:.0f}%)"
            )

    print("Formatting results...")

    # Create a DataFrame from the feature presence dictionary
    combined_features = DataFrame(feature_presence, index=file_names).T

    # Reset index and rename the columns
    combined_features = combined_features.reset_index().rename(
        columns={"index": "Feature"}
    )
    combined_features = combined_features.sort_values(by="Feature")

    # Remove ".tsv" from the column names
    combined_features.columns = [
        col.replace(".tsv", "") for col in combined_features.columns
    ]

    return combined_features


def find_feature_files(features_dir):
    import subprocess

    print("Searching for features files...")

    # Use find because it's much faster than Python os
    command = f"find {features_dir} -name '*.tsv'"
    result = subprocess.check_output(command, shell=True, text=True)

    features_files = sorted(result.strip().split("\n"))

    print(f"Found {len(features_files)} files")

    return features_files


def main():
    """The main script function"""
    from docopt import docopt

    args = docopt(__doc__)

    files = args["<file>"]
    out_file = args["--out-file"]

    if len(files) == 1:
        files = find_feature_files(files[0])
    output = combine_features(files)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
