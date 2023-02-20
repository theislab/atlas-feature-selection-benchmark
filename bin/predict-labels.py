#!/usr/bin/env python

"""
Predict labels for a query dataset

Usage:
    predict-labels.py --reference=<path> --params=<path> --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --reference=<path>   Path to reference dataset file.
    --params=<path>      Path to TSV file containing hyperparameter settings
    --out-file=<path>    Path to output file.
"""


def predict_labels(reference, query, params, seed=1):
    """
    Predict labels for a query dataset

    Parameters
    ----------
    reference
        AnnData object containing the reference dataset
    query
        AnnData object containing the query dataset
    params
        DataFrame containing hyperparameters

    Returns
    -------
    DataFrame containing predicted labels and probabilities for each label class
    """

    from lightgbm import LGBMClassifier
    from sklearn.preprocessing import normalize
    from pandas import DataFrame
    from numpy import argmax

    print("Creating training dataset...")
    X_train = reference.obsm["X_emb"]
    Y_train = reference.obs["Label"].cat.codes
    labels = reference.obs["Label"].cat.categories.tolist()
    n_labels = len(labels)

    print("Using the following parameters...")
    params_dict = dict(zip(params["Parameter"], params["Value"]))
    for param in [
        "n_estimators",
        "num_leaves",
        "max_depth",
        "min_child_samples",
        "subsample_freq",
        "max_bin",
    ]:
        params_dict[param] = int(params_dict[param])

    for param, value in params_dict.items():
        print(f"{param:<25}: {value}")

    print("Training final classification model...")
    lgbm_classifier = LGBMClassifier(
        boosting_type="gbdt",
        metric="multiclass",
        objective="multi_logloss",
        num_class=n_labels,
        n_jobs=1,
        verbose=0,
        random_state=seed,
        **params_dict,
    )
    lgbm_classifier.fit(X_train, Y_train)

    print("Predicting query labels...")
    X_test = query.obsm["X_emb"]
    scores = lgbm_classifier.predict_proba(X_test)
    probs = normalize(scores, axis=1, norm="l1")
    pred_categories = [argmax(line) for line in probs]
    pred_labels = [labels[cat] for cat in pred_categories]

    print("Formatting results...")
    results = DataFrame(
        {
            "ID": query.obs_names,
            "Unseen": query.obs["Unseen"],
            "Label": query.obs["Label"],
            "PredLabel": pred_labels,
            "MaxProb": [
                prob_values[pred_cat]
                for pred_cat, prob_values in zip(pred_categories, probs)
            ],
        }
    )
    for i in range(len(labels)):
        results["Prob" + labels[i]] = probs[:, i]

    return results


def main():
    """The main script function"""
    from docopt import docopt
    from anndata import read_h5ad
    from pandas import read_csv

    args = docopt(__doc__)

    file = args["<file>"]
    reference_file = args["--reference"]
    params_file = args["--params"]
    out_file = args["--out-file"]

    print(f"Reading query data from '{file}'...")
    query = read_h5ad(file)
    print("Read query data:")
    print(query)
    print(f"Reading reference data from '{file}'...")
    reference = read_h5ad(reference_file)
    print("Read reference data:")
    print(reference)
    print(f"Reading parameters from '{file}'...")
    params = read_csv(params_file, sep="\t")
    print("Read parameters:")
    print(params)
    output = predict_labels(reference, query, params)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
