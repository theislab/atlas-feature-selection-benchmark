#!/usr/bin/env python

"""
Predict labels for a query dataset

Usage:
    predict-labels.py --reference=<path> --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --reference=<path>   Path to reference dataset file.
    --out-file=<path>    Path to output file.
"""


def predict_labels(reference, query, seed=1):
    """
    Predict labels for a query dataset

    Parameters
    ----------
    reference
        AnnData object containing the reference dataset
    query
        AnnData object containing the query dataset
    seed
        Random seed to use

    Returns
    -------
    DataFrame containing predicted labels and probabilities for each label class
    """

    from numpy.random import seed as np_seed
    from sklearn.linear_model import LogisticRegression
    from pandas import DataFrame

    np_seed(seed)

    print("Creating training dataset...")
    X_train = reference.obsm["X_emb"]
    Y_train = reference.obs["Label"].cat.codes
    labels = reference.obs["Label"].cat.categories.tolist()

    print("Training classifier...")
    logreg_classifier = LogisticRegression(
        penalty="none",
        class_weight="balanced",
        random_state=seed,
        solver="lbfgs",
        max_iter=10000,
        multi_class="multinomial",
        verbose=1,
        n_jobs=-1,
    )
    logreg_classifier.fit(X_train, Y_train)

    print("Creating test dataset...")
    X_test = query.obsm["X_emb"]
    print("Predicting labels...")
    Y_pred = logreg_classifier.predict(X_test)
    Y_pred_proba = logreg_classifier.predict_proba(X_test)

    print("Formatting results...")
    results = DataFrame(
        {
            "ID": query.obs_names,
            "Unseen": query.obs["Unseen"],
            "Label": query.obs["Label"],
            "PredLabel": [labels[i] for i in Y_pred],
            "MaxProb": [max(p) for p in Y_pred_proba],
        }
    )
    for i, label in enumerate(labels):
        results["Prob" + label] = Y_pred_proba[:, i]

    return results


def main():
    """The main script function"""
    from docopt import docopt
    from anndata import read_h5ad

    args = docopt(__doc__)

    file = args["<file>"]
    reference_file = args["--reference"]
    out_file = args["--out-file"]

    print(f"Reading query data from '{file}'...")
    query = read_h5ad(file)
    print("Read query data:")
    print(query)
    print(f"Reading reference data from '{file}'...")
    reference = read_h5ad(reference_file)
    reference.obs["Label"] = reference.obs["Label"].cat.remove_unused_categories()
    print("Read reference data:")
    print(reference)
    output = predict_labels(reference, query)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
