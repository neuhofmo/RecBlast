# /usr/bin/env python2

import pickle
from sys import exit
from RecBlastUtils import *

# Description:
# This is part 2 of the RecBlast stand-alone version.
# It starts from the first blast results and ends with running the second blast.


# PARAMETERS #
COUNT = 0


def increment():
    """Increments a counter. A global function which is used throughout the script."""
    global COUNT
    COUNT += 1


def gi_to_fasta(awaiting_second_blast, accession_list, matching_orgs, rbh_dict_org, DEBUG, debug, attempt_no=0):
    """
    The function fetches fasta sequences using Entrez.efetch, parses them and returns a dictionary with the parsed data.
    :param awaiting_second_blast: a dictionary (empty or from a previous iteration of the function)
    :param accession_list: list of accession numbers
    :param matching_orgs: a dictionary with matching organisms
    :param rbh_dict_org: a dictionary with organisms with results
    :param DEBUG: a flag
    :param debug: the debug function
    :param attempt_no:  the current attempt (0 by default).
    :return:
    """
    try:  # connecting to ENTREZ nuc DB

        # handle = Entrez.efetch(db="nucleotide", id=accession_list, retmode="xml")  # changed 31/1/17
        handle = Entrez.efetch(db="protein", id=accession_list, retmode="xml")  # changed 31/1/17
        records = Entrez.read(handle)
        debug("successfully connected")

    except Exception, e:  # DB connection exception
        print "Error connecting to server, trying again..."
        print "Error: {}".format(e)
        debug("Error connecting to server, trying again...\n")

        # sleeping in case it's a temporary database problem
        sleep_period = 180
        print "restarting attempt in {} seconds...".format(sleep_period)
        sleep(sleep_period)

        # counting the number of attempts to connect.
        attempt_no += 1
        if attempt_no >= 10:  # If too many:
            print "Tried connecting to Entrez DB more than 10 times. Check your connection or try again later."
            exit(1)

        # try again (recursive until max)
        return gi_to_fasta(awaiting_second_blast, accession_list, matching_orgs, rbh_dict_org, DEBUG, debug, attempt_no)

    # parsing the fasta sequences we retrieved
    for record in records:
        try:
            increment()  # updating COUNT
            accession_v = record["GBSeq_accession-version"]
            try:  # using the old header format (until Sept. 2016)
                fasta_header = "|".join([">GI", record["GBSeq_other-seqids"][1].split("|")[1],
                                         accession_v, record["GBSeq_definition"], record["GBSeq_source"]])
            except IndexError:  # using new format
                fasta_header = "|".join([">", accession_v, record["GBSeq_definition"], record["GBSeq_source"]])
            fasta_record = "\n".join([fasta_header, record["GBSeq_sequence"]])
            # checking if it's RBH using the rbh_dict
            if rbh_dict_org[matching_orgs[accession_v]] == accession_v:
                rbh = True
            else:
                rbh = False
            # save the parsed results to the dictionary
            awaiting_second_blast[COUNT] = (fasta_record, matching_orgs[accession_v], rbh)
            debug("Added {} to awaiting_second_blast dictionary.".format(COUNT))
        except Exception, e:
            print "A problem happened with retrieving the sequences.\nSkipping, but this is the error: {}".format(e)
    return awaiting_second_blast


def split_blast_files(blast_out_filename, local_dic, first_blast_folder, DEBUG, debug):
    """ Going over the blast output file, and splitting it."""
    with open(blast_out_filename, 'r') as blast_out:
        for line in blast_out:
            # getting the sequence identifier (which is the key for the local_id_dic)
            local_seq_id = strip(split(line, '\t')[-1])
            # debug("local_seq_id is {}".format(local_seq_id))  # DEBUG print
            try:
                # reading the blast file output name
                blast_out_name = local_dic[local_seq_id][1]
                # debug("blast_out_name: {}".format(blast_out_name))  # DEBUG print
                full_blast_out_name = join_folder(first_blast_folder, blast_out_name)
                # adding this line to the original file
                with open(full_blast_out_name, 'a') as specific_blast_out:
                    specific_blast_out.write(line)
                    debug("Wrote the following line to {0}:\n{1}".format(full_blast_out_name, line))
            except Exception, e:
                print "Problem with parsing line, sequence identifier may not match the dictionary. Check it."
                print line
                print str(e)
    return True


def main(local_id_dic, first_blast_folder, second_blast_folder, original_id, e_value_thresh, identity_threshold,
         coverage_threshold, accession_regex, run_folder, blastp_path, target_db, outfmt, max_target_seqs,
         back_e_value_thresh, cpu, org_tax_id, run_all, DEBUG, debug, run_second_blast, input_list=None):
    """
    Main function of part 2. Starts from first blast and ends with the second blast.
    """

    if input_list:
        assert(type(input_list) == list)  # make sure it's a list
        listing = input_list
    else:  # get the input from the blast folder.
        listing = sorted([join_folder(first_blast_folder, x) for x in os.listdir(first_blast_folder)
                          if x.endswith('.taxa_filtered.txt')])

    # manually define the id of the gene you are starting with (default is 0)

    blast_two_output_files = []  # initializing update_match_results list for second blast
    blast_two_gene_id_paths = []  # folder paths for the second blast
    second_blast_for_ids_dict = {}  # a dictionary containing all the parsed genes that will be sent to 2nd blast
    rbh_dict = {}  # a dictionary where key=original_id, value={key=organism: value=accession}

    if run_all:  # at this point I added the unified files
        debug("Splitting BLAST results to individual files.")
        filtered_blast_out_filename = join_folder(first_blast_folder, "all_results.taxa_filtered.txt")
        if split_blast_files(filtered_blast_out_filename, local_id_dic, first_blast_folder, DEBUG, debug):
            debug("Blast output parsed successfully.")  # The files should be split by now and we are converging
        else:
            debug("Parsing did not work. exiting.")
            exit(1)

    # each file represents blast results for one gene.
    # per each file, we inspect every line and check if its quality is high enough to become a
    # candidate for the second blast
    for blast_result_file in listing:
        # a successful candidate will be added to this list
        awaiting_second_blast = {}  # dictionary of putative orthologues for this specific gene ID
        # will be blasted by the end of this loop!
        rbh_dict[original_id] = {}  # defining the rbh dictionary for the gene

        # each file represents a different gene, thus different id
        global COUNT
        COUNT = 0  # initialize COUNT
        print "working on file %s" % blast_result_file
        debug("working on file %s , ID: %s\n" % (blast_result_file, original_id))

        # assigning folder for the blast results
        new_blast_path = join_folder(second_blast_folder, "id_{}".format(original_id))
        if os.path.exists(new_blast_path):
            debug("Folder {} already exists!".format(new_blast_path))
        else:
            os.mkdir(new_blast_path)  # creating new folder for the file
            debug("Created new blast first_blast_folder at {}".format(new_blast_path))
        blast_two_gene_id_paths.append(new_blast_path)

        try:
            with open(blast_result_file, 'r') as f:  # open file for reading
                debug("iterating over lines in %s" % blast_result_file)

                # initialize:
                successful_match_gis = []  # for each gene (=file), add the seq_id of every good match
                matching_organisms = {}  # to match the organism in the successful accession_ids

                for line in f:  # read every line
                    # parse every line and take relevant information
                    line_list = split(line)
                    try:  # if it's the old format: gi|761916618|ref|XP_011407293.1|
                        accession_v = split(line_list[1], '|')[3]
                    except IndexError:  # new format: using regex
                        result = re_search(accession_regex, line_list[2])
                        accession_v = result.group(1)
                    identity = line_list[2]
                    coverage = line_list[3]
                    evalue = line_list[4]
                    organism = " ".join([line_list[5], line_list[6]])

                    debug("seq_id: {0}; identity: {1}; coverage: {2}; evalue: {3}; organism: {4}".format(accession_v,
                                                                                                         identity,
                                                                                                         coverage,
                                                                                                         evalue,
                                                                                                         organism))

                    # quality check. if the candidate is good enough, we keep it
                    if float(evalue) < e_value_thresh \
                            and float(identity) >= identity_threshold \
                            and float(coverage) >= coverage_threshold:
                        successful_match_gis.append(accession_v)  # append good match to successful_match_gis
                        matching_organisms[accession_v] = organism
                        # the idea is to save the organism the first time we see it to identify RBH
                        try:
                            debug("Already seen organism {0} for id {1}".format(rbh_dict[original_id][organism],
                                                                                original_id))
                        except KeyError:
                            rbh_dict[original_id][organism] = accession_v

                print "found %s matches for second blast!" % len(successful_match_gis)

                # If there were no matches - moving on to the next line
                if len(successful_match_gis) == 0:
                    print "No successful matches at all for original id {0} in file {1}".format(original_id,
                                                                                                blast_result_file)
                    original_id += 1  # increment
                    continue  # moving on to next file

                # we can't send long lists of genes to entrez, so if the list is long we send it in chunks
                # per each candidate for second blast we extract FASTA and info about the taxa it belongs to
                elif len(successful_match_gis) > 200:
                    chunks = [successful_match_gis[x:x + 100] for x in xrange(0, len(successful_match_gis), 100)]
                    for chunk in chunks:
                        awaiting_second_blast_partial = gi_to_fasta(awaiting_second_blast, chunk, matching_organisms,
                                                                    rbh_dict[original_id], DEBUG, debug)
                        debug("Done importing chunk.")
                        # merge all the dictionaries after splitting to chunks
                        awaiting_second_blast = merge_two_dicts(awaiting_second_blast, awaiting_second_blast_partial)

                else:  # do it for all the genes without splitting
                    awaiting_second_blast = gi_to_fasta(awaiting_second_blast, successful_match_gis, matching_organisms,
                                                        rbh_dict[original_id], DEBUG, debug)
                    debug("Done converting all gis for id %s" % original_id)
        except IOError:
            debug("No results for file {}".format(blast_result_file))

        # documenting all the genes awaiting second blast, per gene
        pickle.dump(awaiting_second_blast, open(join_folder(run_folder, 'genes_for_inspection_id_%s.p' % original_id),
                                                'wb'))
        # adding to dict to be used later:
        second_blast_for_ids_dict[original_id] = awaiting_second_blast

        print "Preparing fasta files and running blast:"

        files_to_blast = []  # initialize fasta files to perform blast on
        debug("awaiting_second_blast: {}".format(awaiting_second_blast))
        for match_id, fasta_seq_list in awaiting_second_blast.iteritems():
            # each item is a gene we need to create a fasta file for

            # file names:
            job_name = "id_%s_index_%s" % (original_id, match_id)  # the job name (also used in the file name)
            fasta_name = "secondblast_fasta_for_{}.fasta".format(job_name)
            full_fasta_path = join_folder(new_blast_path, fasta_name)
            debug("Going to run {0} in full first_blast_folder: {1}".format(job_name, full_fasta_path))
            blast_output_file = replace(full_fasta_path, ".fasta", "_full.txt")
            filtered_blast_output_file = replace(full_fasta_path, ".fasta", ".taxa_filtered.txt")
            debug("blast_output_file: {}".format(blast_output_file))

            # write fasta to file:
            with open(full_fasta_path, "w") as output_fasta_file:
                output_fasta_file.write(fasta_seq_list[0])
            # add the fasta file to the list of files to run the blast commands on
            files_to_blast.append(full_fasta_path)

            # command line to run:
            command_line = "{0} -query {1} -db {2} -outfmt '{3}' -max_target_seqs {4} -evalue {5} -max_hsps 1 " \
                           "-qcov_hsp_perc {6} -num_threads {7} -out {8}\n" \
                           "grep -v ';' {8} | grep {9} > {10}\nrm {8}\n".format(blastp_path, full_fasta_path,
                                                                                target_db, outfmt, max_target_seqs,
                                                                                back_e_value_thresh, coverage_threshold,
                                                                                cpu, blast_output_file, org_tax_id,
                                                                                filtered_blast_output_file)
            if run_second_blast:
                debug("Running the following line:\n{}".format(command_line))

                # writing the command to file and running the file
                try:                                                    # this try paragraph was added later to handle
                    script_path = write_blast_run_script(command_line, run_folder)  # I/O operations,
                    subprocess.check_call(script_path)                              # delay in read/write operations...
                except subprocess.CalledProcessError:                   # restarting the process (with a little sleep)
                    debug("Had a little problem with running this command: "
                          "{}\nSo we are running it again.".format(command_line))
                    sleep(10)
                    script_path = write_blast_run_script(command_line, run_folder)
                    sleep(20)
                    subprocess.check_call(script_path)

            else:
                debug("Not running blast.\nSkipped running the following line:\n{}".format(command_line))
            # adding the filtered update_match_results file name here:
            blast_two_output_files.append(filtered_blast_output_file)

        print "Ran second BLASTP jobs for ID {}".format(original_id)

        original_id += 1  # increment

    debug("Done with part 2.")
    # dumping second blast and rbh dict as backup
    pickle.dump(second_blast_for_ids_dict, open(join_folder(run_folder, "second_blast_for_ids_dict.p"), 'wb'))
    pickle.dump(rbh_dict, open(join_folder(run_folder, "rbh_dict.p"), 'wb'))
    return second_blast_for_ids_dict, blast_two_output_files, blast_two_gene_id_paths

# Done with part 2
