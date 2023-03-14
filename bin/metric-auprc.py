#!/usr/bin/env python

"""
Evaluate cell label classification using the Area Under the Precision-Recall Curve (AUPRC) score

Usage:
    metric-auprc.py --dataset=<str> --method=<str> --integration=<str> --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --out-file=<path>    Path to output file.
"""


def calculate_auprc(labels):
    """
    Calculate AUPRC score for a set of cell labels

    Parameters
    ----------
    labels
        DataFrame containing real and predicted cell labels

    Returns
    -------
    The AUPRC score
    """

    from sklearn.metrics import average_precision_score

    print("Calculating AUPRC score for each label...")
    label_categories = labels["Label"].unique()
    label_scores = []
    for label in label_categories:
        label_true = (labels["Label"] == label).astype(int)
        label_pred = labels["Prob" + label]
        score = average_precision_score(label_true, label_pred)
        print(f"Label '{label}': {score:.4f}")
        label_scores.append(score)

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
    print("Removing unseen populations...")
    input = input[~input["Unseen"]]
    print("Read data:")
    print(input)
    score = calculate_auprc(input)
    output = format_metric_results(
        dataset, method, integration, "Classification", "AUPRC", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
