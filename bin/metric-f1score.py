#!/usr/bin/env python

"""
Evaluate cell label classification using the F1 score

Usage:
    metric-f1score.py --dataset=<str> --method=<str> --integration=<str> [--average=<str>] --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --average=<str>      The type of averaging performed when there are multiple labels [default: "micro"].
    --out-file=<path>    Path to output file.
"""


def calculate_f1score(labels, average=None):
    """
    Calculate F1 score for a set of cell labels

    Parameters
    ----------
    labels
        DataFrame containing real and predicted cell labels
    average
        The type of averaging performed when there are multiple labels. See https://scikit-learn.org/stable/modules/generated/sklearn.metrics.f1_score.html.
        Can also be set to "rarity" to weight labels according to how uncommon they are.

    Returns
    -------
    The F1 score
    """
    from sklearn.metrics import f1_score

    if average == "rarity":
        from pandas import DataFrame

        print("Calculating label weights...")
        label_stats = DataFrame({"Count": labels["Label"].value_counts()})
        label_stats["Frequency"] = label_stats["Count"] / labels.shape[0]
        inverse_freq_sum = sum(1 / label_stats["Frequency"])
        label_stats["Weight"] = 1 / (label_stats["Frequency"] * inverse_freq_sum)

        print("Calculating label scores...")
        label_stats["Score"] = f1_score(
            labels["Label"], labels["PredLabel"], labels=label_stats.index, average=None
        )
        label_stats["WeightedScore"] = label_stats["Weight"] * label_stats["Score"]

        print(label_stats)

        print("Calulating weighted score...")
        score = sum(label_stats["WeightedScore"])
    else:
        score = f1_score(labels["Label"], labels["PredLabel"], average=average)

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
    average = args["--average"]
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    input = read_csv(file, sep="\t")
    print("Removing unseen populations...")
    input = input[~input["Unseen"]]
    print("Read data:")
    print(input)
    score = calculate_f1score(input, average=average)
    output = format_metric_results(
        dataset, method, integration, "Classification", f"F1-{average}", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
