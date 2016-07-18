#! /usr/bin/env python2

from RecBlastUtils import *
import urllib2
import pickle
from time import strftime
import subprocess

# TODO: edit description
# THIS SCRIPT GENERATES FASTA AND SGE FILES FOR THE GENES WE WANT TO USE FOR THE RECIPROCAL BLAST.
# INPUT: CSV FILE CONTAINING A LIST OF GENES WE WANT TO INSPECT
# OUTPUT: FASTA FILES + SGE FILES, ALL READY TO START RUNNING BLAST


# this dictionary holds information about the genes we set out to check. structure:
# {protein_inner_id (index: 0,1,2..) = [common_id, full_id, uniprot_id]}
id_dic = {}  # initialize


# this functions returns a FASTA sequence according to uniprot id.
def get_uni(uni_id, contact):
    """
    Receives a sequence ID and an email address (optional) and returns the FASTA sequence.
    Connecting to UNIPROT
    :param contact: is the contact details (email address) for using the app. By default it will be the recblast email.
    :param uni_id: uniprot id
    :return:
    """
    url2 = "http://www.uniprot.org/uniprot/"  # uniprot to fasta
    url_uniprot = url2 + uni_id + ".fasta"  # the UNIPROT url
    request = urllib2.Request(url_uniprot)
    request.add_header('User-Agent', 'Python %s' % contact)  #
    response = urllib2.urlopen(request)  # get request
    page = response.read(200000)  # read (up to 200000 lines)
    page2 = page.replace("\n", "|")  # edit response
    return page2


# gene_list_to_csv function
def main(file_path, contact, run_folder, fasta_path, first_blast_folder, fasta_output_folder, blastp_path, db,
         taxa_list_file, outfmt, max_target_seqs, e_value_thresh, coverage_threshold, cpu, DEBUG, debug):

    # defined in advance for efficiency:
    regex = re.compile(r'>.*=\d?\|')
    gene_line_regex = re.compile(r'([A-Za-z0-9]+),(.+),([A-Za-z0-9]+)')

    # initialize a list for the blast update_match_results file paths
    blast_one_output_files = []

    debug("Starting to work on the input csv file:")

    with open(file_path) as f:  # open .csv file
        # initialize index:
        csv_line_index = 0  # csv line
        gene_id_index = 1  # matching genes

        for line in f:
            if csv_line_index > 0:  # skipping the header
                try:
                    # generates a FASTA sequence from each protein in the input CSV
                    gene_line_res = re_search(gene_line_regex, strip(line))  # using regex to search
                    common_id = gene_line_res.group(1)  # gene name
                    full_id = gene_line_res.group(2)  # description
                    uniprot_id = strip(replace(gene_line_res.group(3), "\n", ""))

                    # getting the FASTA sequence from uniprot using uniprot_id
                    fa = get_uni(uniprot_id, contact)
                    # cleaning the FASTA sequence
                    result = re.match(regex, fa)
                    fa = re.sub(regex, '', fa)
                    fa = replace(fa, "\n", "")
                    fa = replace(fa, "|", "")
                    sequence = fa  # before we add the header to the fasta, we want to keep the sequence itself
                    fa = "\n".join([result.group(), fa])  # header and sequence
                    # building a dictionary of the proteins we are going to check:
                    # {protein_inner_id: [fasta, common_id, full_id, uniprot_id]}
                    id_dic[gene_id_index] = [fa, common_id, full_id, uniprot_id, sequence]  # save in a dictionary
                    gene_id_index += 1  # moving to the next gene on the list

                except Exception, e:  # in case uniprot doesn't work, please notify!
                    print "There is a problem with retrieving the sequences from UniProt. " \
                          "Please try again later.\n{}".format(e)
                    # TODO: there also might be a problem with the sequence they entered so check how to verify this
            csv_line_index += 1  # next row

    # creating a file with the information about the genes we are checking in this run.
    # this pickle will be used for reference later
    pickle.dump(id_dic, open(os.path.join(run_folder, 'genes_for_inspection_full.p'), 'wb'))
    debug("Success in updating genes_for_inspection file")

    # generating FASTA files and performing the blast:
    for gene_id_index, valueList in id_dic.iteritems():
        # defining file paths:
        debug("gene_index_id: {}".format(gene_id_index))
        job_name = "fasta-%s-%s" % (gene_id_index, valueList[1])  # job_name: index-common_id
        debug("job_name: {}".format(job_name))
        fasta_filename = os.path.join(fasta_path, "{}.fasta".format(job_name))
        blast_out_filename = "{}_full.txt".format(job_name)  # BLAST update_match_results file
        # the update_match_results file after filtering taxa:
        filtered_blast_out_filename = "{}.taxa_filtered.txt".format(job_name)
        with open(fasta_filename, 'w') as output:
            output.write("{}\n\n".format(valueList[0]))  # write fasta to update_match_results file

        # currently using my path on jekyl for running blast
        blast_output_file = os.path.join(first_blast_folder, blast_out_filename)
        filtered_blast_out_filename = os.path.join(first_blast_folder, filtered_blast_out_filename)

        # copy the fasta file to the fasta_output folder
        fasta_output_filename = os.path.join(fasta_output_folder, os.path.basename(fasta_filename))
        shutil.copy(fasta_filename, fasta_output_filename)  # copy

        # command line to run:
        command_line = "{0} -query {1} -db {2} -outfmt '{3}' -max_target_seqs {4} -evalue {5} -max_hsps 1 " \
                       "-qcov_hsp_perc {6} -num_threads {7} -out {8}\n" \
                       "grep -v ';' {8} | grep -w -f {9} > {10}\nrm {8}\n".format(blastp_path, fasta_filename, db,
                                                                                  outfmt, max_target_seqs,
                                                                                  e_value_thresh, coverage_threshold,
                                                                                  cpu, blast_output_file,
                                                                                  taxa_list_file,
                                                                                  filtered_blast_out_filename)

        debug("Running the following line:\n{}".format(command_line))

        # writing the command to file and running the file
        script_path = write_blast_run_script(command_line)
        subprocess.check_call(script_path)  # TODO: test update_match_results

        # adding the filtered file name here:
        blast_one_output_files.append(filtered_blast_out_filename)

        print "Finished running {0}.".format(job_name)

    print "Prepared and ran first BLAST on all FASTA files."
    # dumping id_dic file for pickle:
    pickle.dump(id_dic, open(os.path.join(run_folder, "id_dic.p"), 'wb'))
    print "Part 1 done at {}".format(strftime('%H:%M:%S'))
    return id_dic, blast_one_output_files


# if __name__ == "__main__":
#     if main():
#         exit(0)
#     else:
#         exit(1)  # if you fail, return 1

# DONE with part 1
