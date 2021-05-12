#! /usr/bin/env python3
# coverage.py

"""
Iterates through a list of previously generated samtools depth command files
to extract ABO gene sequence depth of coverage data, calculates average sample
depth per 100 base intervals, and writes depth information for each sample per
line in a text file.


"""

import os


def main():
    """
    Iterates through files to extract depth of coverage data and writes
    average depth per 100 bases to a text file.
    """
    # path to text file containing sample file names
    path_2_file_list = "/Users/jonathan_stevens/ABO/depth_out.txt"
    # create list of sample file names
    file_list = create_list_of_depth_files(path_2_file_list)
    # create list of sample names
    # create dictionary of sample keys with depth per 100 base list values
    sample_list, depth_dict = process_files_from_list(file_list)
    # generate header for output text file
    header_line = create_file_headers()
    # path to outfile
    out_file = '1000G_100bp_avg_read_depth_of_coverage.txt'
    # print data to file
    print_data_2_file(out_file, header_line, sample_list, depth_dict)


def create_list_of_depth_files(text_file):
    """
    creates a list of file names from a text file
    :param text_file: str, name of text file containing file list
    :return: list, list of file names
    """
    file_list = list()
    with open(text_file, 'r') as tf:
        for line in tf:
            line = line.strip()
            file_list.append(line)
    return file_list


def process_files_from_list(file_list):
    """
    Opens each file in file_list, extracts sample name and depth data, creates
    100 base intervals, and calculates their average values.
    :param file_list: list, list of file names
    :return: sample_list, list of sample names
    :return: depth_dict, dictionary of sample keys with list of depth per 100
    base interval values
    """
    sample_list = list()
    depth_dict = dict()
    for file in file_list:
        # append sample name to sample_list
        sample_name = file.split('_')[0]
        sample_list.append(sample_name)
        # specify path to individual depth file
        basepath = '/Users/jonathan_stevens/ABO/1000G_data/depth/'
        file_path = os.path.join(basepath, file)
        # create file handle for depth file
        file_handle = open(file_path, 'r')
        # extract depth data from file
        sample_depth_list = _get_depth_list(file_handle)
        file_handle.close()
        # split depth data into 100 base interval lists
        depth_intervals_list = _create_depth_intervals_list(sample_depth_list)
        # calculate average value for each interval
        avg_depth_per_interval = \
            _calculate_avg_depth_per_interval(depth_intervals_list)
        # add average depth per interval to dictionary
        depth_dict[sample_name] = avg_depth_per_interval
    return sample_list, depth_dict


def _get_depth_list(file_handle):
    """
    Extracts sample depth of coverage data from column 3 of the text file and
    returns a list of the values
    :param file_handle: file object in read mode containing depth data
    :return: sample_depth_list, a list of depth values from file
    """
    sample_depth_list = list()
    for line in file_handle:
        line = line.split()
        # create list of depth values from column 3 of file
        sample_depth_list.append(int(line[2]))
    return sample_depth_list


def _create_depth_intervals_list(sample_depth_list):
    """
    Splits the depth list into 100 element lists
    :param sample_depth_list: list of depth values
    :return: list of 100 element lists
    """
    depth_intervals_list = list()
    # divide list into 100 element lists
    for i in range(0, len(sample_depth_list), 100):
        depth_intervals_list.append(sample_depth_list[i:i + 100])
    return depth_intervals_list


def _calculate_avg_depth_per_interval(depth_intervals_list):
    """
    Calculates the average depth per interval
    :param depth_intervals_list: list of depths per 100 bases
    :return: avg_depth_per_interval, list of average depths per 100 bases
    """
    avg_depth_per_interval = list()
    # calculate average of each interval
    for interval in depth_intervals_list:
        avg_depth_per_interval.append(round(sum(interval) / len(interval)))
    return avg_depth_per_interval


def create_file_headers():
    """
    Creates header line for output text file
    :return: header_line, list of strings
    """
    header_line = list()
    header_line.append("Sample")
    # generate 100 base interval column headers
    for i in range(133255176, 133385146, 100):
        header_line.append(f'chr9:{i}')
    return header_line


def print_data_2_file(out_file, header_line, sample_list, depth_dict):
    """
    Prints average depth of coverage per 100 bases for all samples to a text
    file
    :param out_file: string, name of out put file
    :param header_line: list, list of header values
    :param sample_list: list, list of sample names
    :param depth_dict: dictionary, sample keys with list of depth values
    :return: None
    """
    with open(out_file, 'w') as fh:
        # write tab delimited header line to file
        header_line_joined = "\t".join(header_line)
        fh.write(f'{header_line_joined}\n')
        # write sample data per line to file
        for sample in sample_list:
            coverage = '\t'.join(str(i) for i in depth_dict[sample])
            sample_line = '\t'.join((sample, coverage))
            fh.write(f'{sample_line}\n')


if __name__ == "__main__":
    main()
