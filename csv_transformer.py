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
    :param gene_list: A list of gene names, uniprot IDs or accession numbers.
    :param taxid: The reference taxon ID
    :param out_file: Output file path.
    :return:
    """
    # update_match_results file
    with open(out_file, "w") as output:
        # generating the csv we need for the analysis
        output.write("{}\n".format(",".join(["gene_id", "gene_name", "uniprot_id"])))  # write title
        mg = mygene.MyGeneInfo()  # mygene module
        genes_data = mg.querymany(gene_list, scopes='symbol,name,reporter,accession,ensemblprotein,retired',
                                  fields='uniprot,name,symbol',
                                  species=taxid)

        gene_count = 0
        for geneDic in genes_data:  # iterating over dicts
            try:  # parsing
                original_gene_id = str(geneDic['symbol'])
                full_gene_name = str(geneDic['name'])
                uniprot_id = str(geneDic['uniprot']['Swiss-Prot'])
                output.write("{}\n".format(",".join([original_gene_id, full_gene_name, uniprot_id])))
            except KeyError:  # retrieving didn't succeed
                print("didn't find value for %s" % geneDic['query'])  # didn't happen so far but still
                print("Trying to get the value manually from uniprot in a patchy way:")
                try:  # if it's a uniprot ID, this part should work.
                    url2 = "http://www.uniprot.org/uniprot/"  # uniprot to fasta
                    url_uniprot = url2 + geneDic['query'] + ".tab"  # the UNIPROT url
                    request = urllib2.Request(url_uniprot)
                    response = urllib2.urlopen(request)  # get request
                    page = response.read(20000)  # read (up to 200000 lines)
                    this_gene_data = split(split(page, '\n')[1], '\t')
                    original_gene_id = this_gene_data[1]
                    full_gene_name = this_gene_data[3]
                    uniprot_id = this_gene_data[0]
                    output.write("{}\n".format(",".join([original_gene_id, full_gene_name, uniprot_id])))
                except urllib2.HTTPError:
                    print("didn't find value for %s in uniprot too." % geneDic['query'])

            gene_count += 1

    if os.stat(out_file).st_size > 30:  # minimum file size because of the header
        # this means that the CSV creation worked.
        return True
    else:
        print("Could not convert gene names to csv file. Couldn't locate the genes. Exiting.")
        exit(1)
        return False  # we actually don't need to return False because of the exit. but still.


def gene_file_to_csv(infile, tax_id, outfile=None):
    """
    Receive a file with a list of genes, as well as an optional output file name.
    Returns the filename of the update_match_results csv file.
    :param infile: A file containing a list of genes.
    :param tax_id: Tax ID of the reference genome
    :param outfile: Output file path
    :return:
    """
    # if the outfile is not provided we will generate our own
    if not outfile:
        outfile = "{}.genes.csv".format(infile)
    input_accession_list = gene_list_from_file(infile)
    if gene_list_to_csv(input_accession_list, tax_id, outfile):
        return outfile
    else:
        return ""  # return an empty file name

if __name__ == "__main__":
    from sys import argv
    # list of gene names or accession, input by user
    try:
        input_file = argv[1]
        tax_id = argv[2]
    except IndexError:
        print "Please provide a file containing a list of gene names or accession numbers, and a tax_id (number)."
        exit(1)
    try:
        output_file = argv[2]
    except IndexError:
        output_file = None
    output_file = gene_file_to_csv(input_file, tax_id, output_file)
    print "Wrote CSV to {}".format(output_file)


