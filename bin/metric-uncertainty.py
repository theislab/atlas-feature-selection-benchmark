#!/usr/bin/env python

"""
Evaluate unseen populations using the classification uncertainty metric

Usage:
    metric-uncertainty.py --dataset=<str> --method=<str> --integration=<str> --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --out-file=<path>    Path to output file.
"""


def calculate_uncertainty(labels):
    """
    Calculate the uncertainty score for unseen labels

    Parameters
    ----------
    labels
        DataFrame containing predicted cell labels and probabilities for unseen
        populations

    Returns
    -------
    The unseen uncertainty score
    """

    print("Calculating label scores...")
    unseen_labels = labels["Label"].unique()
    label_scores = []
    for  label in unseen_labels:
        label_probs = labels[labels["Label"] == label]["MaxProb"]
        label_score = 1 - label_probs.mean()
        print(f"Label '{label}': {label_score:.4f}")
        label_scores.append(label_score)

    print("Calculating final score...")
    score = sum(label_scores) / len(label_scores)

    return score


def main():
    """The main script function"""
    from docopt import docopt
    from pandas import read_csv
    from functions.metrics import format_metric_results

    args = docopt(__doc__)
    file = args["<file>"]
    dataset = args["--dataset"]
    method = args["--method"]
    integration = args["--integration"]
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    input = read_csv(file, sep="\t")
    print("Selecting unseen populations...")
    input = input[input["Unseen"]]
    print("Read data:")
    print(input)
    score = calculate_uncertainty(input)
    output = format_metric_results(
        dataset, method, integration, "Unseen", "Uncertainty", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
