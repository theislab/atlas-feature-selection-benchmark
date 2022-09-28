#!/usr/bin/env python

"""
Evaluate cell label classification using the jaccard index metric

Usage:
<<<<<<<< HEAD:bin/metric-f1score.py
    metric-f1score.py --dataset=<str> --method=<str> --integration=<str> [--average=<str>] --out-file=<path> [options] <file>
	
========
    metric-JaccardIndex.py --dataset=<str> --method=<str> --integration=<str> [--average=<str>] --out-file=<path> [options] <file>

>>>>>>>> main:bin/metric-jaccardIndex.py
Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
<<<<<<<< HEAD:bin/metric-f1score.py
    --average=<str>      The parameter is required for multiclass/multilabel targets, if its multiclass, default setting is 'weighted', otherwise is "binary" [default: "None"].
========
    --average=<str>      The parameter is required for multiclass/multilabel targets, if its multiclass, default setting is 'samples', otherwise is "binary" [default: "None"].
>>>>>>>> main:bin/metric-jaccardIndex.py
    --out-file=<path>    Path to output file.
"""


<<<<<<<< HEAD:bin/metric-f1score.py
def calculate_f1score(labels,average=None):
========
def calculate_JaccardIndex(labels,average=None):
>>>>>>>> main:bin/metric-jaccardIndex.py
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
<<<<<<<< HEAD:bin/metric-f1score.py
    from sklearn.metrics import f1_score
			
    score = f1_score(list(labels["Label"]), list(labels["PredLabel"]),average=average)
========
    from sklearn.metrics import jaccard_score

    score = jaccard_score(list(labels["Label"]), list(labels["PredLabel"]),average=average)
>>>>>>>> main:bin/metric-jaccardIndex.py

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
    score = calculate_f1score(input,average=average)
    output = format_metric_results(
<<<<<<<< HEAD:bin/metric-f1score.py
        dataset, method, integration, "Classification", "F1 Score", score
========
        dataset, method, integration, "Classification", f"JaccardIndex{average}", score
>>>>>>>> main:bin/metric-jaccardIndex.py
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
