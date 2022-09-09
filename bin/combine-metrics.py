#!/usr/bin/env python

"""
Combine a set of metrics output files

Usage:
    combine-metrics.py --out-file=<path> [options] <file>...

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
"""


def combine_metrics(files):
    """
    Combine a set of metrics output files

    Parameters
    ----------
    files
        List of paths to metric output files

    Returns
    -------
    DataFrame containing the combined metrics
    """

    from pandas import read_csv, concat

    print(f"Reading metrics from {len(files)} files...")
    metrics_list = (read_csv(file, sep="\t") for file in files)
    print("Combining metrics...")
    metrics = concat(metrics_list, ignore_index=True)
    metrics.sort_values(
        ["Dataset", "Method", "Integration", "Type", "Metric"], inplace=True
    )

    return metrics


def main():
    """The main script function"""
    from docopt import docopt

    args = docopt(__doc__)

    files = args["<file>"]
    out_file = args["--out-file"]

    output = combine_metrics(files)
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
