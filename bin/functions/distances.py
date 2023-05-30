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
