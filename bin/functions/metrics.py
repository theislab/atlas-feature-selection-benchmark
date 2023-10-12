def format_metric_results(dataset, method, integration, metric_type, metric, value):
    """
    Format metric results

    Parameters
    ----------
    dataset
        The name of the dataset the metric has been calculated for
    method
        The name of the method the metric has been calculated for
    integration
        The name of the integration the metric has been calculated for
    metric_type
        The type of the metric, either 'IntegrationBatch', 'IntegrationBio',
        'Classification', 'Mapping' or 'Unseen'
    metric
        The name of the metric that has been calculated
    value
        The value of the calculated metric

    Returns
    -------
    DataFrame containing the formatted results
    """

    from pandas import DataFrame

    if not metric_type in [
        "IntegrationBatch",
        "IntegrationBio",
        "Classification",
        "Mapping",
        "Unseen",
    ]:
        raise ValueError(
            "'metric_type' must be one of 'IntegrationBatch', 'IntegrationBio', 'Classification', 'Mapping' or 'Unseen'"
        )

    if metric ==  "NA":
        from warnings import warn
        warn(f"Storing missing score for metric '{metric}'")
    else:
        if value < 0 or value > 1:
            raise ValueError("'score' must be between 0 and 1")

    print("Formatting metric results...")
    results = DataFrame(
        [
            {
                "Dataset": dataset,
                "Method": method,
                "Integration": integration,
                "Type": metric_type,
                "Metric": metric,
                "Value": value,
            }
        ]
    )

    return results

