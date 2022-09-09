#!/usr/bin/env python

"""
Download the scIB Pancreas dataset (DOI: 10.1038/s41592-021-01336-8)

Usage:
    template.py --out-file=<path> [options]

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
"""


def get_pancreas():
    """
    Get the scIB Pancreas dataset

    Returns
    -------
    AnnData containing the Pancreas dataset
    """

    from scanpy import read
    from tempfile import TemporaryDirectory
    from os.path import join

    url = "https://figshare.com/ndownloader/files/24539828"
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

    output = get_pancreas()
    print("Read dataset:")
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.write_h5ad(out_file)
    print("Done!")


if __name__ == "__main__":
    main()
