#!/usr/bin/env python

"""
Combine a set of metrics output files. If only one file is given it is assumed
that is a directory to search for TSV metric output files.

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

    total_files = len(files)
    progress_step = total_files // 10

    print(f"Reading metrics from {total_files} files...")

    metrics = DataFrame()
    missing_values = DataFrame(columns=["Integration", "Metric", "Missing"])
    missing_files = []

    for i, file_path in enumerate(files, 1):
        try:
            metric_data = read_csv(file_path, sep="\t")

            # Check if the metric value is missing
            is_value_missing = metric_data["Value"].isna().any()

            if is_value_missing:
                integration = metric_data["Integration"].iloc[0]
                metric = metric_data["Metric"].iloc[0]
                existing_row = missing_values[
                    (missing_values["Integration"] == integration)
                    & (missing_values["Metric"] == metric)
                ]

                if not existing_row.empty:
                    missing_values.loc[existing_row.index, "Missing"] += 1
                else:
                    new_row = DataFrame(
                        {
                            "Integration": [integration],
                            "Metric": [metric],
                            "Missing": [1],
                        }
                    )
                    missing_values = concat(
                        [missing_values, new_row], ignore_index=True
                    )

                missing_files.append(file_path)

            # Add to combined metrics
            metrics = concat([metrics, metric_data], ignore_index=True)

            if i % progress_step == 0 or i == total_files:
                print(
                    f"Progress: {i}/{total_files} files processed ({i/total_files * 100:.0f}%)"
                )

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
        integration_summary = missing_values.groupby("Integration")["Missing"].sum()
        print(integration_summary)

        print("\nSummary by Metric:")
        metric_summary = missing_values.groupby("Metric")["Missing"].sum()
        print(metric_summary)

    return (metrics, missing_values, missing_files)


def find_metrics_files(metrics_dir):
    import subprocess

    print("Searching for metrics files...")

    # Use find because it's much faster than Python os
    command = f"find {metrics_dir} -name '*.tsv' ! -name 'all-metrics.tsv' ! -name 'missing-summary.tsv'"
    result = subprocess.check_output(command, shell=True, text=True)

    metrics_files = result.strip().split("\n")

    print(f"Found {len(metrics_files)} files")

    return metrics_files


def main():
    """The main script function"""
    from docopt import docopt

    args = docopt(__doc__)

    files = args["<file>"]
    out_file = args["--out-file"]
    missing_values_file = args["--missing-values"]
    missing_files_file = args["--missing-files"]

    if len(files) == 1:
        files = find_metrics_files(files[0])
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
