#!/usr/bin/env python

"""
Evaluate mapping using the reconstruction error metric

Usage:
    metric-reconstruction.py --dataset=<str> --method=<str> --integration=<str> --exprs=<path> --out-file=<path> <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --exprs=<path>       Path to the query expression H5AD file.
    --out-file=<path>    Path to output file.
"""

# Code in this file is based on https://github.com/MarioniLab/oor_benchmark/blob/8a7a0d41b90d052e5fe2f47db38244041b3896ba/src/oor_benchmark/methods/scArches_mappingQC.py
# available under the MIT license:
#
# MIT License
#
# Copyright (c) 2022, Emma Dann
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import scvi


def calculate_reconstruction_error(model):
    """
    Calculate the reconstruction error score for a mapped dataset.

    Parameters
    ----------
    model
        A trained scvi-tools integration model

    Returns
    -------
    The [0, 1] score for query-reference mixing.
    """

    from scanpy.preprocessing import normalize_total, log1p
    from scipy.sparse import issparse
    from sklearn.metrics.pairwise import cosine_distances
    import numpy as np

    print("Calculating ground-truth expression profiles...")
    adata_truth = model.adata.copy()
    normalize_total(adata_truth, target_sum=10000)
    log1p(adata_truth)
    X_true = adata_truth.X.copy()
    del adata_truth

    print("Sampling from the posterior...")
    # Get 50 samples to get an average reconstruction
    n_samples = 50
    post_samples = np.zeros((model.adata.shape[0], model.adata.shape[1]))
    for i in range(n_samples):
        print(f"Sample {i + 1} of {n_samples}")
        # Get the posterior sample
        post = model.posterior_predictive_sample(n_samples=1)
        while np.any(post.sum(1) == 0):
            print("Posterior has cells with zero counts, resampling...")
            post = model.posterior_predictive_sample(n_samples=1)
        # Normalise the sample
        post = ((post.T / post.sum(1)).T) * 10000
        post = np.log1p(post)
        # Add it to the sum over samples
        post_samples = post_samples + post

    print("Calculating mean reconstructions...")
    X_pred = post_samples / n_samples

    print("Calculating cosine distances...")
    if issparse(X_true):
        X_true = X_true.toarray()
    if issparse(X_pred):
        X_pred = X_pred.toarray()
    cosine_dists = cosine_distances(X_true, X_pred)

    print("Calculating final reconstruction score...")
    # Take the diagonal as we are only interested in distances between the same cells
    cell_scores = np.diag(cosine_dists)
    # Take the mean over cells
    # Divide by 2 to get a score in [0, 1] (range of cosine distances is [0, 2])
    # Substract from 1 so that higher is better
    score = 1 - (np.mean(cell_scores) / 2)

    return score


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad
    from os.path import dirname
    from functions.metrics import format_metric_results

    args = docopt(__doc__)

    query_file = args["<file>"]
    dataset = args["--dataset"]
    method = args["--method"]
    integration = args["--integration"]
    exprs_file = args["--exprs"]
    out_file = args["--out-file"]

    print(f"Reading expression data from '{exprs_file}'...")
    exprs = read_h5ad(exprs_file)
    print("Read data:")
    print(exprs)
    model_dir = dirname(query_file)
    print(f"Reading integration model from '{model_dir}'...")
    model_adata = read_h5ad(query_file)
    model_adata.X = exprs[:, model_adata.var_names].X.copy()
    del exprs
    if "scVI" in integration:
        model = scvi.model.SCVI.load(model_dir, adata=model_adata)
    elif "scANVI" in integration:
        model_adata.obs["ReferenceLabel"] = "Unknown"
        model = scvi.model.SCANVI.load(model_dir, adata=model_adata)
    else:
        raise ValueError(f"Unknown integration method: '{integration}'")
    print("Read model:")
    print(model)
    score = calculate_reconstruction_error(model)
    output = format_metric_results(
        dataset, method, integration, "Mapping", "Reconstruction", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()