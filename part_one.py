#! /usr/bin/env python2

from RecBlastUtils import *
import urllib2
import pickle
from time import strftime
import subprocess

# Description:
# The script performs part 1 of the RecBlast program. Starting from a list of genes and taxa and then moving on to
# saving the sequences and running Blast on them.

# this dictionary holds information about the genes we set out to check. structure:
# {protein_inner_id (index: 0,1,2..) = [common_id, full_id, uniprot_id]}
id_dic = {}  # initialize
local_id_dic = {}


# this functions returns a FASTA sequence according to uniprot id.
def get_uni(uni_id, contact):
    """
    Receives a sequence ID and an email address (optional) and returns the FASTA sequence.
    Connecting to UNIPROT
    :param contact: is the contact details (email address) for using the app. By default it will be the recblast email.
    :param uni_id: uniprot id
    :return:
    """
    url = "http://www.uniprot.org/uniprot/"  # uniprot to fasta
    url_uniprot = url + uni_id + ".fasta"  # the UNIPROT url
    request = urllib2.Request(url_uniprot)
    request.add_header('User-Agent', 'Python %s' % contact)  #
    response = urllib2.urlopen(request)  # get request
    page = response.read(200000)  # read (up to 200000 lines)
    page = replace(page, "\n", "|")
    return page


def main(file_path, contact, run_folder, fasta_path, first_blast_folder, fasta_output_folder, blastp_path, db,
         taxa_list_file, outfmt, max_target_seqs, e_value_thresh, coverage_threshold, cpu, run_all, DEBUG, debug,
         run_first_blast):
    """
    Main function of part_one. Performs most
    :param run_first_blast: if True, runs the first blast process (default)
    :param file_path: The gene list file
    :param contact: email address (RecBlast)
    :param run_folder: path to run in
    :param fasta_path: path of fasta files
    :param first_blast_folder: path of first blast files
    :param fasta_output_folder: path to the output fasta files
    :param blastp_path: path to blastp
    :param db: DB location
    :param taxa_list_file: path to the taxa list file
    :param outfmt: blast out format. constant.
    :param max_target_seqs: Blast parameter
    :param e_value_thresh: Blast parameter
    :param coverage_threshold: Blast parameter
    :param cpu: number of CPU to use in BLAST
    :param DEBUG: DEBUG parameter
    :param debug: debug function
    :param run_all: run all sequences as one file (boolean, default False)
    :return:
    """
    # defined in advance for efficiency:
    regex = re.compile(r'>.*=\d?\|')
    # gene_line_regex = re.compile(r'([A-Za-z0-9]+),(.+),([A-Za-z0-9]+)$')  # Probably an earlier version
    gene_line_regex = re.compile(r'([A-Za-z0-9_]+),(.+),([A-Za-z0-9_]+)')

    # initialize a list for the blast output file paths
    blast_one_output_files = []

    debug("Starting to work on the input csv file:")

    with open(file_path) as f:  # open .csv file
        # initialize index:
        csv_line_index = 0  # csv line
        gene_id_index = 1  # matching genes

        for line in f:
            if csv_line_index > 0 and line != '\n':  # skipping the header and empty lines
                try:
                    # generates a FASTA sequence from each protein in the input CSV
                    gene_line_res = re_search(gene_line_regex, strip(line))  # using regex to search
                    common_id = gene_line_res.group(1)  # gene name
                    full_id = gene_line_res.group(2)  # description
                    uniprot_id = strip(replace(gene_line_res.group(3), "\n", ""))

                    # getting the FASTA sequence from uniprot using uniprot_id
                    fa = get_uni(uniprot_id, contact)
                    # cleaning the FASTA sequence
                    result = re_match(regex, fa)
                    fa = re_sub(regex, '', fa)
                    fa = replace(fa, "\n", "")
                    fa = replace(fa, "|", "")
                    grouped_res = result.group()
                    local_seq_id = split(grouped_res, ' ')[0][1:]  # this is the seq_id used by blast
                    debug("local_seq_id is {}".format(local_seq_id))
                    sequence = fa  # before we add the header to the fasta, we want to keep the sequence itself
                    fa = "\n".join([grouped_res, fa])  # header and sequence
                    # building a dictionary of the proteins we are going to check:
                    # {protein_inner_id: [fasta, common_id, full_id, uniprot_id]}
                    id_dic[gene_id_index] = [fa, common_id, full_id, uniprot_id, sequence, local_seq_id]  # added
                    # save in a dictionary
                    gene_id_index += 1  # moving to the next gene on the list

                except Exception, e:  # in case uniprot doesn't work, please notify!
                    print "There is a problem with retrieving the sequences from UniProt, in line:\n{0}" \
                          "Please try again later.\n{1}".format(line, e)
            csv_line_index += 1  # next row

    # creating a file with the information about the genes we are checking in this run.
    # this pickle will be used for reference later
    pickle.dump(id_dic, open(join_folder(run_folder, 'genes_for_inspection_full.p'), 'wb'))
    debug("Success in updating genes_for_inspection file")

    # This part is for running the sequences together and not individually
    if run_all:
        all_fasta_filename = join_folder(fasta_path, "all_fasta.fasta")
        all_blast_output_file = join_folder(first_blast_folder, "all_results.txt")
        filtered_all_blast_out_filename = join_folder(first_blast_folder, "all_results.taxa_filtered.txt")
        # checking if the file already exists from a previous run:
        if exists_not_empty(all_fasta_filename):
            debug("Fasta file {} already exists, deleting and starting a new one.".format(all_fasta_filename))
            os.remove(all_fasta_filename)

    # generating FASTA files and performing the blast:
    for gene_id_index, valueList in id_dic.iteritems():
        # defining file paths:
        debug("gene_index_id: {}".format(gene_id_index))
        job_name = "fasta-%s-%s" % (gene_id_index, valueList[1])  # job_name: index-common_id
        debug("job_name: {}".format(job_name))
        fasta_filename = join_folder(fasta_path, "{}.fasta".format(job_name))
        blast_out_filename = "{}_full.txt".format(job_name)  # BLAST update_match_results file
        # the update_match_results file after filtering taxa:
        filtered_blast_out_filename = "{}.taxa_filtered.txt".format(job_name)
        # debug(valueList[5])  # That was a debug print.
        local_id_dic[valueList[5]] = (fasta_filename, filtered_blast_out_filename)
        with open(fasta_filename, 'w') as output:
            output.write("{}\n\n".format(valueList[0]))  # write fasta to output file

        # Adding writing to the unified
        if run_all:
            with open(all_fasta_filename, 'a') as output:
                output.write("{}\n\n".format(valueList[0]))  # write fasta to output file

        blast_output_file = join_folder(first_blast_folder, blast_out_filename)
        filtered_blast_out_filename = join_folder(first_blast_folder, filtered_blast_out_filename)

        # copy the fasta file to the fasta_output folder
        fasta_output_name = replace(fasta_filename, '.fasta', '')
        fasta_output_filename_rbh = join_folder(fasta_output_folder,
                                                os.path.basename(fasta_output_name) + '_RBH.fasta')
        fasta_output_filename_ns = join_folder(fasta_output_folder,
                                               os.path.basename(fasta_output_name) + '_non-strict.fasta')
        fasta_output_filename_strict = join_folder(fasta_output_folder,
                                                   os.path.basename(fasta_output_name) + '_strict.fasta')
        # 3 fasta output files:
        shutil.copy(fasta_filename, fasta_output_filename_rbh)  # copy
        shutil.copy(fasta_filename, fasta_output_filename_ns)  # copy
        shutil.copy(fasta_filename, fasta_output_filename_strict)  # copy

        if not run_all:
            # command line to run:
            command_line = "{0} -query {1} -db {2} -outfmt '{3}' -max_target_seqs {4} -evalue {5} -max_hsps 1 " \
                           "-qcov_hsp_perc {6} -num_threads {7} -out {8}\n" \
                           "grep -v ';' {8} | grep -w -f {9} > {10}\nrm {8}\n".format(blastp_path, fasta_filename, db,
                                                                                      outfmt, max_target_seqs,
                                                                                      e_value_thresh,
                                                                                      coverage_threshold, cpu,
                                                                                      blast_output_file, taxa_list_file,
                                                                                      filtered_blast_out_filename)
            if run_first_blast:
                debug("Running the following line:\n{}".format(command_line))

                # writing the command to file and running the file
                try:                                                    # this try paragraph was added later to handle
                    script_path = write_blast_run_script(command_line, run_folder)  # I/O operations,
                    subprocess.check_call(script_path)                              # delay in read/write operations
                except subprocess.CalledProcessError:                   # restarting the process
                    debug("Had a little problem with running this command: "  # (with a short sleep period)
                          "{}\nSo we are running it again.".format(command_line))
                    sleep(10)
                    script_path = write_blast_run_script(command_line, run_folder)
                    sleep(20)
                    subprocess.check_call(script_path)

                print "Finished running {0}.".format(job_name)
            else:
                debug("Not running blast.\nSkipped running the following line:\n{}".format(command_line))

        # adding the filtered file name here:
        blast_one_output_files.append(filtered_blast_out_filename)  # adding even if we didn't run blast

    # Running on the
    if run_all:
        # command line to run:
        command_line = "{0} -query {1} -db {2} -outfmt '{3}' -max_target_seqs {4} -evalue {5} -max_hsps 1 " \
                       "-qcov_hsp_perc {6} -num_threads {7} -out {8}\n" \
                       "grep -v ';' {8} | grep -w -f {9} > {10}\nrm {8}\n".format(blastp_path, all_fasta_filename, db,
                                                                                  outfmt, max_target_seqs,
                                                                                  e_value_thresh, coverage_threshold,
                                                                                  cpu, all_blast_output_file,
                                                                                  taxa_list_file,
                                                                                  filtered_all_blast_out_filename)

        debug("Running the following line:\n{}".format(command_line))

        # writing the command to file and running the file
        try:  # this try paragraph was added later to handle
            script_path = write_blast_run_script(command_line, run_folder)  # I/O operations, delay in read/write
            subprocess.check_call(script_path)                              # operations, etc.
        except subprocess.CalledProcessError:  # restarting the process (with a little sleep period)
            debug("Had a little problem with running this command: "
                  "{}\nSo we are running it again.".format(command_line))
            sleep(10)
            script_path = write_blast_run_script(command_line, run_folder)
            sleep(20)
            subprocess.check_call(script_path)

        print "Finished running all sequences."

    print "Prepared and ran first BLAST on all FASTA files."
    # dumping id_dic file for pickle:
    pickle.dump(id_dic, open(join_folder(run_folder, "id_dic.p"), 'wb'))
    print "Part 1 done at {}".format(strftime('%H:%M:%S'))
    return id_dic, blast_one_output_files, local_id_dic

# DONE with part 1
