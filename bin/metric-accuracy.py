#!/usr/bin/env python

"""
Evaluate cell label classification using the accuracy metric

Usage:
    metric-accuracy.py --dataset=<str> --method=<str> --integration=<str> --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --out-file=<path>    Path to output file.
"""


def calculate_accuracy(labels):
    """
    Calculate classification accuracy for a set of cell labels

    Parameters
    ----------
    labels
        DataFrame containing real and predicted cell labels

    Returns
    -------
    The accuracy score
    """
    from sklearn.metrics import accuracy_score

    score = accuracy_score(labels["Label"], labels["PredLabel"])

    return score


def main():
    """The main script function"""
    from docopt import docopt
    from pandas import read_csv
    from _functions import format_metric_results

    args = docopt(__doc__)

    file = args["<file>"]
    dataset = args["--dataset"]
    method = args["--method"]
    integration = args["--integration"]
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    input = read_csv(file, sep="\t")
    print("Read data:")
    print(input)
    score = calculate_accuracy(input)
    output = format_metric_results(
        dataset, method, integration, "Classification", "Accuracy", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()