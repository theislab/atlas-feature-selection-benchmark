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
    from scipy.spatial.distance import cosine as cosine_dist
    from tqdm import tqdm
    from math import ceil
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
    batch_size = 10000
    n_cells = model.adata.n_obs
    n_batches = ceil(n_cells / batch_size)
    cell_scores = np.zeros(0)

    for batch in range(n_batches):
        start = batch * batch_size
        end = min(start + batch_size, n_cells)
        print(f"Sampling cells {start + 1} to {end} (batch {batch + 1} of {n_batches})")
        batch_indices = list(range(start, end))
        batch_n = len(batch_indices)
        post_samples = np.zeros((batch_n, model.adata.shape[1]))

        for _ in tqdm(range(n_samples), desc="Sampling", leave=False):

            post = model.posterior_predictive_sample(
                indices=batch_indices, n_samples=1
            )

            # Normalise the sample
            post = ((post.T / post.sum(1)).T) * 10000
            post[~ np.isfinite(post)] = 0 # Reset all zero cells to 0
            post = np.log1p(post)
            # Add it to the sum over samples
            post_samples = post_samples + post

        print("Calculating mean reconstructions...")
        X_pred = post_samples / n_samples

        print("Calculating cosine distances...")
        batch_scores = np.zeros(batch_n)
        X_true_batch = X_true[batch_indices, :]
        # Calculate the cosine distance between each cell's true and predicted expression profile
        for cell in tqdm(range(batch_n), desc="Calculating distances", leave=False):
            true_cell = X_true_batch[cell, :]
            pred_cell = X_pred[cell, :]
            if issparse(true_cell):
                true_cell = true_cell.toarray().ravel()
            if issparse(pred_cell):
                pred_cell = pred_cell.toarray().ravel()
            batch_scores[cell] = cosine_dist(true_cell, pred_cell)

        cell_scores = np.concatenate((cell_scores, batch_scores))

    print("Calculating final reconstruction score...")
    # Take the mean over cells
    # Divide by 2 to get a score in [0, 1] (range of cosine distances is [0, 2])
    # Substract from 1 so that higher is better
    score = 1 - (np.mean(cell_scores) / 2)

    return score


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad
    from scanpy.preprocessing import filter_cells
    from os.path import dirname
    from functions.metrics import format_metric_results

    args = docopt(__doc__)

    query_file = args["<file>"]
    dataset = args["--dataset"]
    method = args["--method"]
    integration = args["--integration"]
    exprs_file = args["--exprs"]
    out_file = args["--out-file"]

    if "Symphony" in integration:
        from warnings import warn

        warn(
            "Reconstruction metric can not be computed for Symphony integrations. Returning 'NA' score."
        )
        score = "NA"
        output = format_metric_results(
            dataset, method, integration, "Mapping", "Reconstruction", score
        )
        print(output)
        print(f"Writing output to '{out_file}'...")
        output.to_csv(out_file, sep="\t", index=False)
        print("Done!")
        return

    print(f"Reading expression data from '{exprs_file}'...")
    exprs = read_h5ad(exprs_file)
    print("Read data:")
    print(exprs)
    model_dir = dirname(query_file)
    print(f"Reading integration model from '{model_dir}'...")
    model_adata = read_h5ad(query_file)
    model_adata.X = exprs[:, model_adata.var_names].X.copy()
    del exprs
    n_cells_pre = model_adata.n_obs
    filter_cells(model_adata, min_genes=5)
    filter_cells(model_adata, min_counts=10)
    if not model_adata.n_obs == n_cells_pre:
        n_cells_removed = n_cells_pre - model_adata.n_obs
        print(f"\nRemoved {n_cells_removed} cells with less than 5 genes or 10 counts")
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
