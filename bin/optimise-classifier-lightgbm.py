#!/usr/bin/env python

"""
Optimise parameters for a label classifier

Usage:
    optimise-classifier.py --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
"""


def optimise_lightgbm(input, n_folds=10, seed=1):
    """
    Optimise parameters for a LightGBM model using Bayesian Optimization

    Based on https://www.kaggle.com/code/lucamassaron/scikit-optimize-for-lightgbm/notebook

    Parameters
    ----------

    input
        AnnData containing the training data
    n_folds
        The number of cross-validation folds to use
    seed
        Random seed to use

    Returns
    -------
    Dictionary with optimized parameters
    """

    from lightgbm import LGBMClassifier
    from pandas import DataFrame
    from sklearn.metrics import log_loss, make_scorer
    from sklearn.model_selection import StratifiedKFold
    from skopt import BayesSearchCV
    from skopt.callbacks import DeltaYStopper, VerboseCallback
    from skopt.space import Real, Integer

    print("Creating training dataset...")
    X_train = input.obsm["X_emb"]
    Y_train = input.obs["Label"].cat.codes
    labels = sorted(input.obs["Label"].cat.categories.tolist())
    n_labels = len(labels)

    logloss_scorer = make_scorer(
        log_loss,
        greater_is_better=False,
        needs_proba=True,
        labels=sorted(Y_train.unique()),
    )

    skf_splitter = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)

    lgbm_classifier = LGBMClassifier(
        boosting_type="gbdt",
        metric="multiclass",
        objective="multi_logloss",
        num_class=n_labels,
        n_jobs=1,
        verbose=-1,
        random_state=seed,
    )

    search_spaces = {
        "learning_rate": Real(0.001, 1.0, "log-uniform"),  # Boosting learning rate
        "n_estimators": Integer(30, 2000),  # Number of boosted trees to fit
        "num_leaves": Integer(4, 200),  # Maximum tree leaves for base learners
        "max_depth": Integer(
            5, 80
        ),  # Maximum tree depth for base learners, <=0 means no limit
        "min_child_samples": Integer(20, 256),  # Minimal number of data in one leaf
        "max_bin": Integer(
            10, 500
        ),  # Max number of bins that feature values will be bucketed
        "subsample": Real(
            0.1, 1.0, "uniform"
        ),  # Subsample ratio of the training instance
        "subsample_freq": Integer(0, 10),  # Frequency of subsample, <=0 means no enable
        "colsample_bytree": Real(
            0.01, 1.0, "uniform"
        ),  # Subsample ratio of columns when constructing each tree
        "min_child_weight": Real(
            0.0001, 10.0, "uniform"
        ),  # Minimum sum of instance weight (hessian) needed in a child (leaf)
        "reg_lambda": Real(1e-9, 100.0, "log-uniform"),  # L2 regularization
        "reg_alpha": Real(1e-9, 100.0, "log-uniform"),  # L1 regularization
    }

    optimiser = BayesSearchCV(
        estimator=lgbm_classifier,
        search_spaces=search_spaces,
        scoring=logloss_scorer,
        cv=skf_splitter,
        n_iter=100,  # Maximum number of trials
        n_points=1,  # Number of hyperparameter sets evaluated at the same time
        n_jobs=-1,  # Number of jobs
        return_train_score=False,
        refit=False,
        random_state=seed,
        verbose=0,
    )

    delta_stopper = DeltaYStopper(delta=0.01)
    verbosity = VerboseCallback(n_total=1)

    print(f"Exploring {optimiser.total_iterations} parameter sets...")
    optimiser.fit(X_train, Y_train, callback=[delta_stopper, verbosity])

    cv_results = DataFrame(optimiser.cv_results_)
    best_score = optimiser.best_score_
    best_score_std = cv_results.iloc[optimiser.best_index_].std_test_score
    best_params = optimiser.best_params_

    print(f"{len(optimiser.cv_results_['params'])} parameter sets checked")
    print(f"Best score: {-best_score}\u00B1{best_score_std:.2f}")

    print("Best parameters:")
    for param, value in best_params.items():
        print(f"{param:<25}: {value}")

    output = DataFrame({"Parameter": best_params.keys(), "Value": best_params.values()})

    return output


def main():
    """The main script function"""
    from docopt import docopt
    from anndata import read_h5ad

    args = docopt(__doc__)

    file = args["<file>"]
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    input.obs["Label"] = input.obs["Label"].cat.remove_unused_categories()
    print("Read data:")
    print(input)
    output = optimise_lightgbm(input)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
