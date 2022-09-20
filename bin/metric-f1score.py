#!/usr/bin/env python

"""
Evaluate cell label classification using the f1 score metric

Usage:
    metric-f1score.py --dataset=<str> --method=<str> --integration=<str> [--average=<str>] --out-file=<path> [options] <file>
	
Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
	--average=<str>      The parameter is required for multiclass/multilabel targets, if its multiclass, default setting is 'weighted', otherwise is "binary" [default: "None"].
    --out-file=<path>    Path to output file.
"""


def calculate_f1score(labels,average="None"):
    """
    Calculate classification accuracy for a set of cell labels
    Parameters
    ----------
    labels
        DataFrame containing real and predicted cell labels
	average
		the parameters required for multiclass targets
    Returns
    -------
    The F1 score
    """
    from sklearn.metrics import f1_score
	
    if average == "None":
        len_mark = len(set(list(labels["Label"])))
        if len_mark == 2:
            average = "binary"
        else:
            average = "weighted"
			
    score = f1_score(list(labels["Label"]), list(labels["PredLabel"]),average=average)

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
    average = str(average)

    print(f"Reading data from '{file}'...")
    input = read_csv(file, sep="\t")
    print("Read data:")
    print(input)
    score = calculate_f1score(input,average=average)
    output = format_metric_results(
        dataset, method, integration, "Classification", "F1 Score", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()