#!/usr/bin/env python

"""
Select features using simple statistics

Usage:
    method-statistic.py --out-file=<path> [options] <file>

Options:
    -h --help               Show this screen.
    -o --out-file=<path>    Path to output file.
    -n --n-features=<int>   Number of features to select [default: 2000].
    -f --statistic=<str>    Statistic to use for feature selection. One of: 'mean', 'variance' [default: mean].
"""


def select_features_statistic(adata, n_features, statistic):
    """
    Select features using simple statistics

    Parameters
    ----------
    adata
        AnnData object
    n_features
        Number of features to select
    statistic
        Statistic to use ('mean', 'variance')

    Returns
    ----------
    DataFrame containing the selected features
    """

    from scanpy.preprocessing import normalize_total, log1p

    if adata.n_vars < n_features:
        import warnings

        warnings.warn(
            "Number of features to select is greater than the number present, setting n_features to adata.n_vars"
        )
        n_features = adata.n_vars

    print("Normalising expression...")
    normalize_total(adata, target_sum=1e4)
    log1p(adata)

    if statistic == "mean":
        from scanpy.preprocessing import calculate_qc_metrics

        print("Calculating feature means...")
        _, feature_stats = calculate_qc_metrics(adata, expr_type="logcounts", percent_top=None)
        feature_stats["Feature"] = adata.var_names
        feature_stats = feature_stats.sort_values(
            by="mean_logcounts", ascending=False
        )

    if statistic == "variance":
        from numpy import var
        from pandas import DataFrame

        print("Calculate feature variances...")
        variances = var(adata.X, axis=0)
        feature_stats = DataFrame(adata.var_names, columns=["Feature"])
        feature_stats["Variance"] = variances
        feature_stats = feature_stats.sort_values(
            by="Variance", ascending=False
        )

    print(f"Selecting top {n_features} features...")
    selected_features = feature_stats.head(n=n_features)

    return selected_features


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad

    args = docopt(__doc__)

    file = args["<file>"]
    n_features = int(args["--n-features"])
    statistic = args["--statistic"]
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    input = read_h5ad(file)
    print("Read data:")
    print(input)
    output = select_features_statistic(input, n_features, statistic)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
