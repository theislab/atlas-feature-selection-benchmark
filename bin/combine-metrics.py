#!/usr/bin/env python

"""
Combine a set of metrics output files

Usage:
    combine-metrics.py --out-file=<path> --missing-values=<path> --missing-files=<path> [options] <file>...

Options:
    -h --help                Show this screen.
    --out-file=<path>        Path to combined metrics output file.
    --missing-values=<path>  Path to missing values output file.
    --missing-files=<path>   Path to missing values files output file.
"""


def combine_metrics(files):
    """
    Combine a set of metrics output files

    Parameters
    ----------
    files
        List of paths to metric output files

    Returns
    -------
    Tuple containing:
        - DataFrame containing the combined metrics
        - DataFrame containing counts of missing metric values
        - List of files with missing metric values
    """

    from pandas import read_csv, concat, DataFrame

    print(f"Reading metrics from {len(files)} files...")

    metrics = DataFrame()
    missing_values = DataFrame(columns=["Integration", "Metric", "Missing"])
    missing_files = []

    for file_path in files:
        try:
            metric_data = read_csv(file_path, sep="\t")

            # Check if the metric value is missing
            is_value_missing = metric_data["Value"].isna().any()
            
            if is_value_missing:
                integration = metric_data["Integration"].iloc[0]
                metric = metric_data["Metric"].iloc[0]
                existing_row = missing_values[
                    (missing_values["Integration"] == integration) &
                    (missing_values["Metric"] == metric)
                ]
                
                if not existing_row.empty:
                    missing_values.loc[existing_row.index, "Missing"] += 1
                else:
                    new_row = DataFrame(
                        {
                            "Integration": [integration],
                            "Metric": [metric],
                            "Missing": [1]
                        }
                    )
                    missing_values = concat([missing_values, new_row], ignore_index=True)
                
                missing_files.append(file_path)

            # Add to combined metrics
            metrics = concat([metrics, metric_data], ignore_index=True)
        except Exception as e:
            # Raise an error if reading fails
            raise ValueError(f"Error reading file {file_path}: {str(e)}")

    metrics.sort_values(
        ["Dataset", "Method", "Integration", "Type", "Metric"], inplace=True
    )

    print(metrics)

    n_missing = len(missing_files)
    print(f"\n{n_missing} files with missing metric values")

    if n_missing > 0:
        print(missing_values)
        
        print("\nSummary by Integration:")
        integration_summary = missing_values.groupby('Integration')['Missing'].sum()
        print(integration_summary)

        print("\nSummary by Metric:")
        metric_summary = missing_values.groupby('Metric')['Missing'].sum()
        print(metric_summary)

    return (metrics, missing_values, missing_files)


def main():
    """The main script function"""
    from docopt import docopt

    args = docopt(__doc__)

    files = args["<file>"]
    out_file = args["--out-file"]
    missing_values_file = args["--missing-values"]
    missing_files_file = args["--missing-files"]

    output, missing_values, missing_files = combine_metrics(files)
    print(f"Writing output to '{out_file}'...")
    output.to_csv(out_file, sep="\t", index=False)
    print(f"Writing missing values to '{missing_values_file}'...")
    missing_values.to_csv(missing_values_file, sep="\t", index=False)
    print(f"Writing missing values files to '{missing_files_file}'...")
    with open(missing_files_file, "w") as outfile:
        for filepath in missing_files:
            outfile.write(filepath + "\n")
    print("Done!")


if __name__ == "__main__":
    main()
