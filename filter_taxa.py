#! /usr/bin/env python2

import gzip
import shutil
import os


def filter_taxa(blast_file, taxa_list):
    """
    receives a blast file and a list of taxa IDs.
    Returns a filtered update_match_results file path.
    :param blast_file:
    :param taxa_list:
    :return:
    """

    output_file = "{}.taxa_filtered.txt".format(blast_file[:blast_file.rfind('.')])

    with open(output_file, 'w') as outfile:  # open update_match_results file
        # read blast input file
        with open(blast_file, 'r') as blast_infile:
            for line in blast_infile:  # go over lines
                tax_id = line.split()[0]  # take first column
                if tax_id in taxa_list:  # if it's in the list
                    outfile.write(line)  # write to update_match_results file

    print "Result written to {}".format(output_file)
    return output_file


def filter_all_taxa(blast_file_list, taxa_file=None, taxa_list=None):
    """

    :param taxa_list:
    :param blast_file_list:
    :param taxa_file:
    :return:
    """
    # reading taxa file list
    if taxa_file:
        with open(taxa_file, 'r') as taxa_infile:
            taxa_list = [line.strip() for line in taxa_infile.readlines()]

    output_file_list = []
    for bfile in blast_file_list:  # iterating over update_match_results file list
        filtered_file = filter_taxa(bfile, taxa_list)
        # compressing the old file
        with open(bfile, 'rb') as f_in, gzip.open(bfile + ".gzip", 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        output_file_list.append(filtered_file)
        os.remove(bfile)
        print("Gzipped and removed file {}".format(bfile))
    return output_file_list  # returning list of filtered files

