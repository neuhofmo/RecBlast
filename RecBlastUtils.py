#! /usr/bin/env python2

# A set of tools, functions, aliases and more used in RecBlast.

import os
import tarfile
from time import strftime, sleep
import re
import subprocess
from Bio import Entrez
import shutil


# debugging function
def debug_s(debug_string, to_debug):
    """
    Receives a string and prints it, with a timestamp.
    :param debug_string: a string to print
    :param to_debug: boolean flag: True means print, False - ignore.
    :return:
    """
    if to_debug:
        print "DEBUG {0}: {1}".format(strftime('%H:%M:%S'), debug_string)


def create_folder_if_needed(path):
    """
    Receives a path and creates a folder when needed (if it doesn't already exist).
    """
    if os.path.exists(path):
        print "{} dir exists".format(path)
    else:
        print "{} dir does not exist. Creating dir.".format(path)
        os.mkdir(path)


def targz_list(archive_name, file_list):
    """
    Returns True after
    :param archive_name:
    :param file_list:
    :return:
    """
    with tarfile.open(archive_name, "w:gz") as tar:
        for file_name in file_list:
            tar.add(file_name, arcname=os.path.basename(file_name))
            os.remove(file_name)  # removing the sge file
    return True


def cleanup(path, fasta_path, first_blast_folder, second_blast_folder):
    """
    Performs tar and gzip on sets of files produced by the program.
    Then deletes the files and folders.
    :param path:  # the run_folder
    :param fasta_path:  # the folder containing fasta files
    :param first_blast_folder:  # the folder containing the first blast results
    :param second_blast_folder:  # the folder containing the second blast results
    :return:
    """
    # clean pickles:
    pickles = [os.path.join(path, x) for x in os.listdir(path) if x[-2:] == '.p']
    pickle_archive = os.path.join(path, "pickles.tar.gz")
    if targz_list(pickle_archive, pickles):
        # pickles cleaned
        pass
    # compress all fasta_path
    fasta_files = [os.path.join(fasta_path, x) for x in os.listdir(fasta_path)]
    fasta_archive = os.path.join(path, "fasta_archive.tar.gz")
    if targz_list(fasta_archive, fasta_files):
        shutil.rmtree(fasta_path)
    # compress all first_blast
    first_blast_files = [os.path.join(first_blast_folder, x) for x in os.listdir(first_blast_folder)]
    first_blast_archive = os.path.join(path, "first_blast_archive.tar.gz")
    if targz_list(first_blast_archive, first_blast_files):
        shutil.rmtree(first_blast_folder)
    # compress all second_blast (recursive)
    second_blast_archive = os.path.join(path, "second_blast_archive.tar.gz")
    with tarfile.open(second_blast_archive, "w:gz") as tar:
        tar.add(second_blast_folder, arcname=os.path.basename(second_blast_folder))
        shutil.rmtree(second_blast_folder)
    return True


def write_blast_run_script(command_line):
    """Writing a blast run script, and giving it run permissions."""
    script_path = "/tmp/blastp_run.sh"  # default script location
    with open(script_path, 'w') as script:
        script.write("#! /bin/tcsh\n")
        script.write("# The script is designed to run the following blastp command from RecBlast\n")
        script.write(command_line)
        # run permissions for the script:
    os.chmod(script_path, 0751)
    return script_path


def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z


def is_number(s):
    """
    The script determines if a string is a number or a text.
    Returns True if it's a number.
    """
    try:
        int(s)
        return True
    except ValueError:
        return False


def blastdb_exit():
    """Exiting if we can't find the $BLASTDB on the local machine"""
    print("$BLASTDB was not found! Please set the blast DB path to the right location.")
    print("Make sure blast+ is installed correctly.")
    exit(1)


def exists_not_empty(path):
    """Receives a file path and checks if it exists and not empty."""
    if os.path.exists(path) and os.stat(path).st_size > 0:
        return True
    else:
        return False


def subset_db(tax_id, gi_file_path, db_path, big_db, run_anyway, DEBUG, debug, attempt_no=0):  # TODO call it, doc it
    """
    Subsets a big blast database into a smaller one based on tax_id.
    The function connects to entrez and retrieves gi identifiers of sequences with the same tax_id.

    :param tax_id: The tax_id (string)
    :param gi_file_path: file path of the gi_list file we are creating
    :param db_path: the new db path
    :param big_db: we are about to subset
    :param run_anyway: run on NR if unable to subset
    :param attempt_no: counter for the attempts in connecting to Entrez (attempts to connect up to 10 times).
    :param DEBUG: A boolean flag: True for debug prints, False for quiet run.
    :param debug: A function call to provide debug prints.
    :return:
    """

    # connecting to ENTREZ protein DB
    try:
        handle = Entrez.esearch(db="protein", term="txid{}[ORGN]".format(tax_id), retmode="xml", retmax=10000000)
        record = Entrez.read(handle)

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
        return subset_db(tax_id, gi_file_path, db_path, big_db, run_anyway, DEBUG, debug, attempt_no)

    assert int(record["Count"]) == len(record["IdList"]), "Did not fetch all sequences!"  # make sure we got it all...
    # writing a gi list file
    with open(gi_file_path, 'w') as gi_file:
        gi_file.write("\n".join(record["IdList"]) + "\n")

    # the new target database path
    create_folder_if_needed(os.path.join(db_path, tax_id))
    target_db = os.path.join(db_path, tax_id, "db")
    aliastool_command = ["blastdb_aliastool", "-gilist", gi_file_path, "-db", big_db, "-dbtype", "prot", "-out",
                         target_db]
    try:
        subprocess.check_call(aliastool_command)
        print("Created DB subset from nr protein for {}".format(tax_id))
        return target_db
    except subprocess.CalledProcessError:
        print("Problem with creating DB for tax_id {} from nr.".format(tax_id))
        if run_anyway:
            print("Running with the heavy nr option. Do some stretches, it might be a long run.")
            return big_db
        print("Aborting.\n"
              "If you want to run the program anyway against the entire nr "
              "(which is significantly slower than the default run, please use the --run_even_if_no_db_found flag.")
        exit(1)


# for efficiency
strip = str.strip
split = str.split
replace = str.replace
re_search = re.search
upper = str.upper
