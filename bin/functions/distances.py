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

    import numpy

    inverse_covariances = {}
    for group in set(groups):
        is_group = [g == group for g in groups]
        group_coords = coords[is_group, :]
        covariance = numpy.cov(group_coords, rowvar=False)

        # Calculate the inverse covariance
        try:
            inv_covariance = numpy.linalg.inv(covariance)
        except numpy.linalg.LinAlgError:
            # If this fails (i.e. the matrix is singular) calculate the pseudo-inverse instead.
            print(
                f"Warning: Unable to calculate the inverse covariance matrix for group '{group}'. Calculating the pseudo-inverse instead."
            )
            inv_covariance = numpy.linalg.pinv(covariance)

        # Check if the inverse covariance matrix has negative values in the diagonal
        # If so, add a correction value to the diagonal
        # This should only happen in rare cases where coordinates are very close to each other
        if numpy.any(numpy.diag(inv_covariance) < 0):
            print(
                f"Warning: Negative values in the diagonal of the inverse covariance matrix for group '{group}'. Adding correction."
            )
            correction = -numpy.min(numpy.diag(inv_covariance)) * 1.1
            print(f"Correction value: {correction}:")
            numpy.fill_diagonal(inv_covariance, numpy.diag(inv_covariance) + correction)

        inverse_covariances[group] = inv_covariance

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
