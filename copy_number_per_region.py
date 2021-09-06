#! /usr/bin/env python3
# copy_number_per_region.py

"""
Iterates through a list of files to extract ABO gene sequence depth of coverage
data, calculates gene copy number per gene region, and writes copy number
information for each sample per line to a tab delimited text file.

To run:
python3 copy_number_per_region.py

"""

import os


def main():
    # path to text file containing sample file names
    path_2_file_list = "/Users/jonathan_stevens/ABO/depth_out.txt"
    # create list of sample file names
    file_list = create_list_of_depth_files(path_2_file_list)
    #
    sample_list, sample_cn_dict = extract_region_depths_from_files(file_list)
    out_file = "copy_number_per_region.txt"
    print_data_2_file(out_file, sample_list, sample_cn_dict)


def create_list_of_depth_files(text_file):
    """
    creates a list of file names from a text file
    :param text_file: str, name of text file containing file list
    :return: file_list, list of file names
    """
    file_list = list()
    with open(text_file, 'r') as tf:
        for line in tf:
            line = line.strip()
            file_list.append(line)
    return file_list


def extract_region_depths_from_files(file_list):
    sample_list = list()
    sample_cn_dict = dict()

    for file in file_list:
        # capture sample name from file name
        sample_name = file.split('_')[0]
        sample_list.append(sample_name)
        # path to depth files
        basepath = '/Users/jonathan_stevens/ABO/1000G_data/depth/'
        # individual path to file
        file_path = os.path.join(basepath, file)
        # create file handle for depth file
        file_handle = open(file_path, 'r')

        # initialize gene region lists
        exon7_utr3_depths = list()
        intron6_depths = list()
        exon6_depths = list()
        intron5_depths = list()
        exon5_depths = list()
        intron4_depths = list()
        exon4_depths = list()
        intron3_depths = list()
        exon3_depths = list()
        intron2_depths = list()
        exon2_depths = list()
        intron1_depths = list()
        exon1_utr5_depths = list()
        baseline_depths = list()

        for line in file_handle:
            line = line.strip().split()
            # extract depth data per gene region
            if int(line[1]) in range(133255176, 133256357):
                exon7_utr3_depths.append(int(line[2]))
            elif int(line[1]) in range(133256357, 133257409):
                intron6_depths.append(int(line[2]))
            elif int(line[1]) in range(133257409, 133257543):
                exon6_depths.append(int(line[2]))
            elif int(line[1]) in range(133257543, 133258097):
                intron5_depths.append(int(line[2]))
            elif int(line[1]) in range(133258097, 133258133):
                exon5_depths.append(int(line[2]))
            elif int(line[1]) in range(133258133, 133259819):
                intron4_depths.append(int(line[2]))
            elif int(line[1]) in range(133259819, 133259867):
                exon4_depths.append(int(line[2]))
            elif int(line[1]) in range(133259867, 133261318):
                intron3_depths.append(int(line[2]))
            elif int(line[1]) in range(133261318, 133261375):
                exon3_depths.append(int(line[2]))
            elif int(line[1]) in range(133261375, 133262099):
                intron2_depths.append(int(line[2]))
            elif int(line[1]) in range(133262099, 133262169):
                exon2_depths.append(int(line[2]))
            elif int(line[1]) in range(133262169, 133275162):
                intron1_depths.append(int(line[2]))
            elif int(line[1]) in range(133275162, 133275215):
                exon1_utr5_depths.append(int(line[2]))
            # 5000 base region used for copy number calculations
            elif int(line[1]) in range(133279500, 133284501):
                baseline_depths.append(int(line[2]))
            else:
                pass
        # calculate baseline average to be used in calculations
        baseline_depth_avg = _average(baseline_depths)
        baseline_depth = round(baseline_depth_avg, 2)

        # copy number calculations per regions
        try:
            exon1_utr5 = round(_average(exon1_utr5_depths) /
                               (baseline_depth_avg / 2), 2)
            intron1 = round(_average(intron1_depths) /
                            (baseline_depth_avg / 2), 2)
            exon2 = round(_average(exon2_depths) /
                          (baseline_depth_avg / 2), 2)
            intron2 = round(_average(intron2_depths) /
                            (baseline_depth_avg / 2), 2)
            exon3 = round(_average(exon3_depths) /
                          (baseline_depth_avg / 2), 2)
            intron3 = round(_average(intron3_depths) /
                            (baseline_depth_avg / 2), 2)
            exon4 = round(_average(exon4_depths) /
                          (baseline_depth_avg / 2), 2)
            intron4 = round(_average(intron4_depths) /
                            (baseline_depth_avg / 2), 2)
            exon5 = round(_average(exon5_depths) /
                          (baseline_depth_avg / 2), 2)
            intron5 = round(_average(intron5_depths) /
                            (baseline_depth_avg / 2), 2)
            exon6 = round(_average(exon6_depths) /
                          (baseline_depth_avg / 2), 2)
            intron6 = round(_average(intron6_depths) /
                            (baseline_depth_avg / 2), 2)
            exon7_utr3 = round(_average(exon7_utr3_depths) /
                               (baseline_depth_avg / 2), 2)
            # add sample copy number per region to dictionary
            sample_cn_dict[sample_name] = [baseline_depth, exon1_utr5, intron1,
                                           exon2, intron2, exon3, intron3,
                                           exon4, intron4, exon5, intron5,
                                           exon6, intron6, exon7_utr3]
        # check for empty depth files (zero for all values)
        except ZeroDivisionError:
            print(f'{sample_name}: ZeroDivisionError, check depth file')
        file_handle.close()
    return sample_list, sample_cn_dict


def _average(depth_list):
    return sum(depth_list) / len(depth_list)


def print_data_2_file(out_file, sample_list, sample_cn_dict):
    with open(out_file, 'w') as fh:
        # write tab delimited header line to file
        header_line = "\t".join(["Sample", "Baseline_depth", "5'UTR_Exon1",
                                 "Intron1", "Exon2", "Intron2", "Exon3",
                                 "Intron3", "Exon4", "Intron4", "Exon5",
                                 "Intron5", "Exon6", "Intron6", "Exon7_3'UTR"])
        fh.write(f'{header_line}\n')
        # write sample data per line to file
        for sample in sample_list:
            try:
                copy_number = '\t'.join(
                    str(i) for i in sample_cn_dict[sample])
                sample_line = '\t'.join((sample, copy_number))
                fh.write(f'{sample_line}\n')
            except KeyError:
                print(f'{sample}: KeyError, could not write to file')


if __name__ == "__main__":
    main()
