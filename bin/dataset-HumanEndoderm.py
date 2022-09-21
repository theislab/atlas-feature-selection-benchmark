#!/usr/bin/env python

"""
Download the human endoderm atlas (DOI: 10.1016/j.cell.2021.04.028)

Usage:
    dataset-HumanEndoderm.py --out-file=<path> [options]

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
"""


def get_humanEndoderm():
    """
    Get the human endoderm atlas

    Returns
    -------
    AnnData containing the dataset
    """

    from scanpy import read
    from tempfile import TemporaryDirectory
    from os.path import join

    url = "https://md-datasets-cache-zipfiles-prod.s3.eu-west-1.amazonaws.com/x53tts3zfr-1.zip"
    call = join("wget ", url)
    os.system(call)

    url = "https://md-datasets-cache-zipfiles-prod.s3.eu-west-1.amazonaws.com/x53tts3zfr-1.zip"
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

    output = get_humanEndoderm()
    print("Read dataset:")
    print(output)
    print(f"Writing output to '{out_file}'...")
    output.write_h5ad(out_file)
    print("Done!")


if __name__ == "__main__":
    main()
