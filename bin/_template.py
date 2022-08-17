#!/usr/bin/env python

"""
A Python script template

Usage:
    template.py --out-file=<path> [options] <file>

Options:
    -h --help            Show this screen.
    --out-file=<path>    Path to output file.
"""


def run_analysis(input):
    """
    A function that performs analysis

    Parameters
    ----------
    input
        The function input

    Returns
    -------
    The function output
    """

    output = input

    return output


def main():
    """The main script function"""
    from docopt import docopt

    args = docopt(__doc__)

    file = args["<file>"]
    out_file = args["--out-file"]

    print(f"Reading data from '{file}'...")
    # input = read_data(file)
    print("Read data:")
    # print(input)
    output = run_analysis(input)
    print(f"Writing output to '{out_file}'...")
    # write_output(output, out_file)
    print("Done!")


if __name__ == "__main__":
    main()
