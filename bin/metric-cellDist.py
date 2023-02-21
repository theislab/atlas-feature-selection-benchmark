#!/usr/bin/env python

"""
Calculate the cell mapping distance metric

Usage:
    metric-cellDist.py --dataset=<str> --method=<str> --integration=<str> --reference=<path> --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --dataset=<str>      Name of the dataset to calculate the metric for.
    --method=<str>       Name of the method to calculate the metric for.
    --integration=<str>  Name of the integration to calculate the metric for.
    --reference=<path>   Path to the reference H5AD file.
    --out-file=<path>    Path to output file.
"""


def calculate_cell_distance(reference, query, upper_bound=30.0):
    """
    Calculate the cell Mahalanobis mapping distance metric

    Parameters
    ----------
    reference
        H5AD object containing the integrated reference dataset
    query
        H5AD object containing the mapped query dataset
    upper_bound
        The upper bound for cell distances, distances greater than this will be
        clipped to this value to calculate the final score

    Returns
    -------
    The cell mapping distance score
    """

    from scipy.spatial.distance import mahalanobis
    from numpy import mean

    reference_coords = reference.obsm["X_emb"]
    reference_labels = reference.obs["Label"].tolist()

    print("Calculating inverse covariance matrices...")
    inverse_covariances = get_inverse_covariances(reference_coords, reference_labels)
    print("Calculating centroid positions...")
    centroids = get_centroids(reference_coords, reference_labels)

    print("Calculating cell Mahalonobis distances...")
    distances = []
    for idx in range(query.n_obs):
        query_label = query.obs["Label"][idx]
        # Skip labels that are not in the reference
        if query_label not in inverse_covariances.keys():
            continue

        query_coord = query.obsm["X_emb"][idx, :]

        distance = mahalanobis(query_coord, centroids[query_label], inverse_covariances[query_label])
        distances.append(distance)

    print("Calculating cell mapping distance score...")
    if max(distances) > upper_bound:
        from warnings import warn
        warn(f"Some distances are greater than the upper bound ({upper_bound}), these will be clipped to calculate the final score")
        distances = [min(d, upper_bound) for d in distances]

    # To get the final score:
    # - Calculate the mean distance across cells
    # - Subtract this from `upper_bound` as we want lower distances to be better
    # - Divide by `upper_bound` to get score in the range [0, 1]
    score = (upper_bound - mean(distances)) / upper_bound

    return score


def get_inverse_covariances(coords, groups):
    """
    Calculate the inverse covariance matrices for each group

    Parameters
    ----------
    coords
        Array giving the coordinates of each sample
    group
        List giving the group for each sample

    Returns
    -------
    Dictionary where keys are groups and values are inverse covariance matrices
    """

    from numpy import cov
    from numpy.linalg import inv

    inverse_covariances = {}
    for group in set(groups):
        is_group = [g == group for g in groups]
        group_coords = coords[is_group, :]
        covariance = cov(group_coords, rowvar=False)
        inverse_covariances[group] = inv(covariance)

    return inverse_covariances


def get_centroids(coords, groups):
    """
    Calculate mean centroid positions each group

    Parameters
    ----------
    coords
        Array giving the coordinates of each sample
    group
        List giving the group for each sample

    Returns
    -------
    Dictionary where keys are groups and values are centoid positions
    """

    from numpy import mean

    centroids = {}
    for group in set(groups):
        is_group = [g == group for g in groups]
        group_coords = coords[is_group, :]
        centroids[group] = mean(group_coords, axis=0)

    return centroids


def main():
    """The main script function"""
    from docopt import docopt
    from scanpy import read_h5ad
    from _functions import format_metric_results

    args = docopt(__doc__)

    query_file = args["<file>"]
    dataset = args["--dataset"]
    method = args["--method"]
    integration = args["--integration"]
    reference_file = args["--reference"]
    out_file = args["--out-file"]

    print(f"Reading query data from '{query_file}'...")
    query = read_h5ad(query_file)
    query.obs["Dataset"] = "Query"
    print("Read data:")
    print(query)
    print(f"Reading reference data from '{reference_file}'...")
    reference = read_h5ad(reference_file)
    reference.obs["Dataset"] = "Reference"
    print("Read data:")
    print(reference)
    score = calculate_cell_distance(reference, query)
    output = format_metric_results(
        dataset, method, integration, "Mapping", "cellDist", score
    )
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == "__main__":
    main()
