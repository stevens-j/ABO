#! /usr/bin/env python3
# reverse_complement_R.py

"""
Reads a bcftools query generated tab delimited text file by line, reverse
complements the variant sequences, and writes them to a new tab delimited
text file.

usage: reverse_complement.py [-h] [-i INFILE] [-o OUTFILE]

Enter infile and outfile names for reverse complementing.

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        path to the file to open
  -o OUTFILE, --outfile OUTFILE
                        name of outfile

"""

import argparse
from Bio.Seq import Seq


def main():
    """
    Reverse complement variant sequences in a bcftools query generated file.
    """
    # get command line arguments
    args = get_cli_args()
    infile = args.infile
    outfile = args.outfile

    # open text file for reading
    infile_handle = get_fh(file=infile, mode='r')
    # iterate over lines of text file to create dictionary of variants
    header_line, sample_list, sample_variant_dict = \
        get_sample_variant_dictionary(infile_handle)
    infile_handle.close()
    # reverse complement variants in dictionary
    reverse_complement_sample_variant_dict = \
        reverse_complement(sample_list, sample_variant_dict)
    # create outfile for writing
    outfile_handle = get_fh(file=outfile, mode='w')
    # write reverse complemented variants to file
    write_reverse_complements_2_file(outfile_handle, header_line, sample_list,
                                     reverse_complement_sample_variant_dict)
    outfile_handle.close()


def write_reverse_complements_2_file(outfile_handle,
                                     header_line, sample_list,
                                     reverse_complement_sample_variant_dict):
    """
    writes reverse complemented sample variants to tab delimited text file
    :param outfile_handle: file handle for outfile
    :param header_line: str, tab separated header line
    :param sample_list: list of sample names
    :param reverse_complement_sample_variant_dict: dict, sample name keys,
    lists of variant seq values
    :return: none
    """
    outfile_handle.write(f'{header_line}')
    for sample in sample_list:
        variant_str = "\t".join(reverse_complement_sample_variant_dict[sample])
        sample_line = f'{sample}\t{variant_str}\n'
        outfile_handle.write(sample_line)


def reverse_complement(sample_list, sample_variant_dict):
    """
    reverse complements sequences in sample_variant_dict
    :param sample_list: list of sample names
    :param sample_variant_dict: dict of sample keys with variant seq
    lists
    :return: reverse_complement_sample_variant_dict
    """
    reverse_complement_sample_variant_dict = {}
    for sample in sample_list:
        variant_list = sample_variant_dict[sample]
        rc_variant_list = []
        for variant in variant_list:
            # remove newline character if present
            variant = variant.strip()
            # determine if variant is indel
            slash = "/"
            if slash in variant:
                indel = variant.split("/")
                # remove quote marks
                # reverse complement sequences
                seq1 = str(Seq(indel[0].strip('"')).reverse_complement())
                seq2 = str(Seq(indel[1].strip('"')).reverse_complement())
                # rejoin indel variants
                rc_indel = [seq1, seq2]
                my_seq = "/".join(rc_indel)
                # append reverse complemented variant to list
                rc_variant_list.append(my_seq)
            else:
                # reverse complement sequence
                my_seq = str(Seq(variant).reverse_complement())
                # append reverse complemented variant to list
                rc_variant_list.append(my_seq)
        # add sample rc_variant_list to reverse_complement_sample_dict
        reverse_complement_sample_variant_dict[sample] = rc_variant_list
    return reverse_complement_sample_variant_dict


def get_sample_variant_dictionary(file_handle):
    """
    Splits each line of the tab delimited text file and assigns each sample
    name as a dictionary key and the variant sequences as dictionary values
    :param file_handle: tab delimited bcftools query file containing sample
    names at index [0] of each row.
    :return: header_line, a tab delimited str
    sample_list, a list of sample names
    sample_variant_dict, sample name keys with variant sequence list
    values.
    """
    header_line = ""
    sample_list = []
    sample_variant_dict = {}
    # read file by line
    for line in file_handle:
        line = line.split("\t")
        # identify sample lines
        if "X" in line[0]:
            sample = line[0]
            # split on . to extract sample name, index 2
            sample = sample.split(".")[2]
            # add sample name to list
            sample_list.append(sample)
            # assign sample/variants to dictionary
            sample_variant_dict[sample] = line[1:]
        else:
            # header line
            line[0] = "Sample"
            header_line = "\t".join(line)
    return header_line, sample_list, sample_variant_dict


def get_fh(file=None, mode=None):
    """
    Opens a file for reading or writing and returns a handle
    :param file: str, file name
    :param mode: str, mode for opening file, "r" or "w"
    :return: file handle
    """
    try:
        handle = open(file, mode)
        return handle
    except IOError:
        raise IOError(f'File not found: "{file}"')
    except ValueError:
        raise ValueError(f'Could not open file: {file} for mode "{mode}"')


def get_cli_args():
    """
    Get command line options with argparse
    :return: instance of argparse arguments
    """
    parser = argparse.ArgumentParser(
        description='Enter infile and outfile names for reverse complementing')
    parser.add_argument('-i', '--infile', dest='infile',
                        type=str, help='path to the file to open')
    parser.add_argument('-o', '--outfile', dest='outfile',
                        type=str, help='name of outfile')
    return parser.parse_args()


if __name__ == "__main__":
    main()
