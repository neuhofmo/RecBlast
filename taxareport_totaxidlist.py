#! /usr/bin/env python2

# from types import *

# taxa report is generated by http://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi
# this script parses the file generated by the link above and returns a list of taxa ID's, ready to be used by recBlast

# if the user enters a list of taxa ids - not need to run this script.
# if the user enters a list of taxa names, he should go to the abovementioned link and generate a file
# I hoped we could access this url without the user knowing, but there's not api


OUTPUT_PATH = "C://Users/Efrat/Documents/Eva/taxa_list.txt"
TAXA_REPORT = "C://Users/Efrat/Documents/Eva/tax_report.txt"
TAX_DB = "/data/DB/taxdump/tax_names.txt"  # TODO: distribute this
tax_list = []


# def convert_tax_to_taxid(TAXA_ID_PROVIDED, OUTPUT_PATH):
#     # extract taxa ids from input file
#     with open(TAXA_ID_PROVIDED) as f:
#         with open(OUTPUT_PATH, "wb") as update_match_results:
#             for line in f:
#                 if "tax" not in line:
#                     tax_number = line.split("|")[3].strip()
#                     update_match_results.write(tax_number + "\n")
#                     tax_list.append(tax_number)
#         print "done"
#     return tax_list


def create_tax_dict(tax_db=TAX_DB):
    """

    :param tax_db:
    :return:
    """
    with open(tax_db, 'r') as db:
        tax_list = [x.strip().split('\t') for x in db.readlines()]
        tax_dict = dict([(x[1], x[0]) for x in tax_list])
    return tax_dict


def convert_tax_to_taxid(tax_dict, TAX_LIST_FILE, output_path=OUTPUT_PATH):  # TODO: maybe save the dict as a pickle?

    strip = str.strip()
    with open(TAXA_LIST_FILE) as f:
        with open(output_path, "wb") as output:
            for line in f:
                line = strip(line)
                try:
                    tax_id = tax_dict[line]
                except NameError:
                    print "Taxon {} does not exist in NCBI tax!".format(line)
                output.write(tax_id + "\n")
                tax_list.append(tax_id)
        print "done"
    return tax_list

