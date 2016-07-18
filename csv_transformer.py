#! /usr/bin/env python2

import mygene
import os
import urllib2
from RecBlastUtils import split


def gene_list_from_file(in_file):
    """
    Receives a file name with a list of gene names or accessions and turns them into a list.
    :param in_file: a file name
    :return: list of gene names / accession numbers
    """
    with open(in_file, 'r') as accession_name:
        file_lines = accession_name.readlines()
        accession_list = [x.strip() for x in file_lines]
    return accession_list


def gene_list_to_csv(gene_list, taxid, out_file):
    """
    Receives a list of gene names/accession numbers and an update_match_results file path.
    :param gene_list:
    :param taxid:
    :param out_file:
    :return:
    """

    # update_match_results file
    with open(out_file, "w") as output:
        # generating the csv we need for the analysis
        mg = mygene.MyGeneInfo()
        genes_data = mg.querymany(gene_list, scopes='symbol,namereporter,accession', fields='uniprot,name,symbol',
                                  species=taxid)  # TODO: check this on other species
        for geneDic in genes_data:  # iterating over dicts
            try:
                original_gene_id = geneDic['symbol']
                full_gene_name = geneDic['name']
                uniprot_id = geneDic['uniprot']['Swiss-Prot']
                output.write("{}\n".format(",".join([original_gene_id, full_gene_name, uniprot_id])))
            except KeyError:
                print "didn't find value for %s" % geneDic['query']  # didn't happen so far but still
    if os.stat(out_file).st_size > 0:
        return True
    else:
        print("Could not convert gene names to csv file. Couldn't locate the genes. Exiting.")
        exit(1)


def gene_file_to_csv(infile, tax_id, outfile=None):
    """
    Receive a file with a list of genes, as well as an optional update_match_results file name.
    Returns the filename of the update_match_results csv file.
    :param infile:
    :param tax_id:
    :param outfile:
    :return:
    """
    if not outfile:
        outfile = "{}.genes.csv".format(infile)
    input_accession_list = gene_list_from_file(infile)
    if gene_list_to_csv(input_accession_list, tax_id, outfile):
        return outfile

if __name__ == "__main__":
    from sys import argv
    # list of gene names or accession, input by user
    try:
        input_file = argv[1]
        tax_id = argv[2]
        # we used 'control_genes_accesion.txt'  as input
    except IndexError:
        print "Please provide a file containing a list of gene names or accession numbers, and a tax_id (number)."
        exit(1)
    try:
        output_file = argv[2]
    except IndexError:
        output_file = None
        # output_file = "checkList.csv"  # FOR TEST
    output_file = gene_file_to_csv(input_file, tax_id, output_file)
    print "Wrote CSV to {}".format(output_file)

# 6087
