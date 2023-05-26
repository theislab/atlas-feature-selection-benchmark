#!/usr/bin/env python

"""
Download the fetal liver hematopoiesis dataset

Usage:
    dataset-fetalLiver.py --out-file=<path> [options]

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
"""


def get_fetal_liver():
    """
    Get the fetal liver hematopoiesis dataset

    Returns
    -------
    AnnData containing the fetal liver dataset
    """

    from scanpy import read
    from tempfile import TemporaryDirectory
    from os.path import join

    url = "https://app.cellatlas.io/fetal-liver/dataset/1/download"
    print(f"Reading dataset from '{url}'...")
    temp_dir = TemporaryDirectory()
    adata = read(join(temp_dir.name, "temp.h5ad"), backup_url=url)

    print("Cleaning up temporary directory...")
    temp_dir.cleanup()

    return adata


def main():
    """The main script function"""
    from docopt import docopt

    args = docopt(__doc__)

    out_file = args["--out-file"]

    output = get_fetal_liver()
    print("Read dataset:")
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.write_h5ad(out_file)
    print("Done!")


if __name__ == "__main__":
    main()
