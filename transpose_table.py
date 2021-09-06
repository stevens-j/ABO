#! /usr/bin/env python3
# transpose_table.py

"""
Transposes a tab delimited text file with pandas and writes the new
file to the working directory.

usage: transpose_table.py [-h] [-i INFILE] [-o OUTFILE]

Opens tab delimited text file as pandas data frame, transposes dataframe,
and writes new text file to the working directory

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        path to input file
  -o OUTFILE, --outfile OUTFILE
                        name of output file

"""

import argparse
import pandas


def main():
    # get command line arguments
    args = get_cli_args()
    infile = args.infile
    outfile = args.outfile

    # read tab delimited file as pandas dataframe
    df = pandas.read_csv(infile, sep='\t')
    # transpose dataframe
    df_t = df.T
    # write transposed dataframe to new file
    df_t.to_csv(outfile, header=False, sep='\t')


def get_cli_args():
    """
    Get command line options with argparse
    :return: instance of argparse arguments
    """
    parser = argparse.ArgumentParser(
        description='Opens tab delimited text file as pandas data frame, '
                    'transposes dataframe, and writes new text file to the '
                    'working directory')
    parser.add_argument('-i', '--infile', dest='infile', type=str,
                        help='path to input file')
    parser.add_argument('-o', '--outfile', dest='outfile', type=str,
                        help='name of output file')
    return parser.parse_args()


if __name__ == "__main__":
    main()
