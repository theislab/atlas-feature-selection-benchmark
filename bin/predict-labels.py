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


def predict_labels(reference, query):
    """
    Predict labels for a query dataset

    Parameters
    ----------
    reference
        AnnData object containing the reference dataset
    query
        AnnData object containing the query dataset

    Returns
    -------
    DataFrame containing predicted labels and probabilities for each label class
    """

    from lightgbm import Dataset, train, early_stopping, log_evaluation, reset_parameter
    from sklearn.preprocessing import normalize
    from pandas import DataFrame
    from numpy import argmax

    print("Creating training dataset...")
    X_train = reference.obsm["X_emb"]
    Y_train = reference.obs["Label"].cat.codes
    labels = reference.obs["Label"].cat.categories.tolist()
    n_labels = len(labels)

    print("Optimising model hyperparameters...")
    best_score, opt_params = optimise_lightgbm(X_train, Y_train, n_labels)
    print(f"Best score: {best_score}")
    print("Selected GBDT model with the following parameters:")
    for param, value in opt_params.items():
        print(f"{param:<25}: {value}")

    print("Training final classification model...")
    dataset_params = {
        "max_bin": int(round(opt_params["max_bin"])),
    }
    params = {
        "objective": "multiclass",
        "metric": "multi_logloss",
        "num_class": n_labels,
        "bagging_freq": 10,
        "learning_rate": max(min(opt_params["base_learning_rate"], 1), 0),
        "num_leaves": int(round(opt_params["num_leaves"])),
        "bagging_fraction": max(min(opt_params["bagging_fraction"], 1), 0),
        "max_depth": int(round(opt_params["max_depth"])),
        "max_bin": int(round(opt_params["max_bin"])),
        "min_data_in_leaf": int(
            round(opt_params["min_prop_data_in_leaf"] * X_train.shape[0])
        ),
        "min_sum_hessian_in_leaf": opt_params["min_sum_hessian_in_leaf"],
        "min_gain_to_split": 0.001,
        "verbose": -1,
    }

    train_data = Dataset(X_train, Y_train, params=dataset_params)
    model = train(
        params,
        train_data,
        num_boost_round=1000,
        valid_sets=[train_data],
        callbacks=[
            early_stopping(stopping_rounds=10, verbose=True),
            log_evaluation(period=100),
            reset_parameter(
                learning_rate=generate_learning_rate_decay(
                    max(min(opt_params["base_learning_rate"], 1), 0),
                    max(min(opt_params["learning_rate_decay"], 1), 0),
                )
            ),
        ],
    )

    print("Predicting query labels...")
    X_test = query.obsm["X_emb"]
    scores = model.predict(X_test)
    probs = normalize(scores, axis=1, norm="l1")
    pred_categories = [argmax(line) for line in probs]
    pred_labels = [labels[cat] for cat in pred_categories]

    print("Formatting results...")
    results = DataFrame(
        {
            "ID": query.obs_names,
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


def optimise_lightgbm(
    X_train, Y_train, n_labels, init_round=15, opt_round=25, n_folds=10, seed=1
):
    """
    Optimise parameters for a LightGBM model using Bayesian Optimization

    Based on https://www.kaggle.com/code/somang1418/tuning-hyperparameters-under-10-minutes-lgbm

    Parameters
    ----------
    X_train
        The data to use for training
    Y_train
        The target labels to use for training
    n_labels
        Number of label categories
    init_round
        The number of steps for random exploration of the parameter space
    opt_round
        The number of steps for Bayesian Optimization
    n_folds
        The number of cross-validation folds to use
    seed
        Random seed to use

    Returns
    -------
    Dictionary with optimized parameters
    """

    from lightgbm import Dataset, cv, early_stopping, log_evaluation, reset_parameter
    from bayes_opt import BayesianOptimization
    from numpy import mean
    from pandas import Series

    def crossvalidate_lightgbm(
        base_learning_rate,
        learning_rate_decay,
        num_leaves,
        bagging_fraction,
        max_depth,
        max_bin,
        min_prop_data_in_leaf,
        min_sum_hessian_in_leaf,
    ):

        dataset_params = {
            "max_bin": int(round(max_bin)),
        }

        train_data = Dataset(
            X_train, Y_train, params=dataset_params, free_raw_data=False
        )

        params = {
            "objective": "multiclass",
            "metric": "multi_logloss",
            "num_class": n_labels,
            "bagging_freq": 10,
            "learning_rate": max(min(base_learning_rate, 1), 0),
            "num_leaves": int(round(num_leaves)),
            "bagging_fraction": max(min(bagging_fraction, 1), 0),
            "max_depth": int(round(max_depth)),
            "max_bin": int(round(max_bin)),
            "min_data_in_leaf": int(round(min_prop_data_in_leaf * X_train.shape[0])),
            "min_sum_hessian_in_leaf": min_sum_hessian_in_leaf,
            "min_gain_to_split": 0.001,
            "verbosity": -1,
        }

        cv_result = cv(
            params,
            train_data,
            num_boost_round=1000,
            nfold=n_folds,
            stratified=True,
            callbacks=[
                early_stopping(stopping_rounds=10, verbose=False),
                log_evaluation(period=0),
                reset_parameter(
                    learning_rate=generate_learning_rate_decay(
                        max(min(base_learning_rate, 1), 0),
                        max(min(learning_rate_decay, 1), 0),
                    )
                ),
            ],
            seed=seed,
        )

        # We want to minimise the loss so maximise the negative loss
        return -mean(cv_result["multi_logloss-mean"])

    bayes_optimizer = BayesianOptimization(
        crossvalidate_lightgbm,
        {
            "base_learning_rate": (0.1, 0.6),
            "learning_rate_decay": (0.001, 0.01),
            "num_leaves": (4, 50),
            "bagging_fraction": (0.5, 1),
            "max_depth": (5, 50),
            "max_bin": (10, 511),
            "min_prop_data_in_leaf": (0.001, 0.05),
            "min_sum_hessian_in_leaf": (1e-10, 100),
        },
        random_state=1,
    )

    bayes_optimizer.maximize(init_points=init_round, n_iter=opt_round)

    model_auc = []
    for model in range(len(bayes_optimizer.res)):
        model_auc.append(bayes_optimizer.res[model]["target"])

    # Return best parameters
    return (
        bayes_optimizer.res[Series(model_auc).idxmax()]["target"],
        bayes_optimizer.res[Series(model_auc).idxmax()]["params"],
    )


def generate_learning_rate_decay(base_rate, decay_rate):
    """
    Generate a learning rate decay function

    Parameters
    ----------
    base_rate
        The base learning rate
    decay_rate
        The rate that learning decays

    Returns
    -------
    A function that returns the learning rate for each iteration
    """

    def decay_function(iteration):
        from math import exp

        learning_rate = base_rate * exp(-decay_rate * iteration)
        learning_rate = max(1e-4, learning_rate)

        return learning_rate

    return decay_function


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
    print("Read data:")
    print(query)
    print(f"Reading reference data from '{file}'...")
    reference = read_h5ad(reference_file)
    print("Read data:")
    print(reference)
    output = predict_labels(reference, query)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
