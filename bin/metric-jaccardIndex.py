#!/usr/bin/env python

"""
Evaluate cell label classification using the jaccard index metric

Usage:
    metric-JaccardIndex.py --dataset=<str> --method=<str> --integration=<str> [--average=<str>] --out-file=<path> [options] <file>
	
Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --average=<str>      The parameter is required for multiclass/multilabel targets, if its multiclass, default setting is 'samples', otherwise is "binary" [default: "None"].
    --out-file=<path>    Path to output file.
"""


def calculate_JaccardIndex(labels,average=None):
    """
    Calculate jaccard index for a set of cell labels
    Parameters
    ----------
    labels
        DataFrame containing real and predicted cell labels
	average
        the parameters required for multiclass targets
    Returns
    -------
    The jaccard index
    """
    from sklearn.metrics import jaccard_score
	
    score = jaccard_score(list(labels["Label"]), list(labels["PredLabel"]),average=average)

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
    average = args["--average"]
    out_file = args["--out-file"]
    if average == "None":
        average = None
    
    print(f"Reading data from '{file}'...")
    input = read_csv(file, sep="\t")
    print("Read data:")
    print(input)
    score = calculate_JaccardIndex(input,average=average)
    output = format_metric_results(
        dataset, method, integration, "Classification", "Jaccard Index", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()