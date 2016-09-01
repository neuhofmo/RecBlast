#! /usr/bin/env python2

# A set of tools, functions, aliases and more used in RecBlast.

import os
import tarfile
from time import strftime, sleep
import re
import subprocess
from Bio import Entrez
import shutil
import zipfile


TEMP_FILES_PATH = os.getcwd()


def prepare_files(items, file_name, user_id, files_path=TEMP_FILES_PATH):
    """Receives a list of items and a file to write them to, then writes them to file and returns the file path."""
    full_path = os.path.join(files_path, "_".join([user_id, file_name]))
    # items = list(set(items))  # make the list unique  # unnecessary
    with open(full_path, 'w') as f:
        for item in items:
            f.write(item + "\n")
    return full_path


def file_to_string(file_name):
    with open(file_name, 'r') as f:
        text = f.read()
    # delete original file
    os.remove(file_name)
    return text


def remove_commas(file_name):
    """Replaces commas with newlines in a file."""
    with open(file_name, 'r') as f:
        text = f.read()
        text = replace(text, ',', '\n')
    with open(file_name, 'w') as f:  # now writing
        f.write(text)
    return file_name


def zip_results(fasta_output_path, csv_rbh_output_filename, csv_strict_output_filename, csv_ns_output_filename,
                output_path):
    """
    Receives a folder containing fasta sequences and a csv file, adds them all to zip.
    :param fasta_output_path:
    :param csv_rbh_output_filename:
    :param csv_strict_output_filename:
    :param csv_ns_output_filename:
    :param output_path:
    :return:
    """
    zip_file = os.path.join(output_path, "output.zip")
    fastas = [os.path.join(fasta_output_path, x) for x in os.listdir(fasta_output_path)]
    with zipfile.ZipFile(zip_file, mode='w') as zf:
        # adding all fasta files
        for fasta in fastas:
            zf.write(fasta, os.path.basename(fasta))
        # zf.write(csv_file_path, os.path.basename(csv_file_path))  # add csv file
        # add csv files
        zf.write(csv_rbh_output_filename, os.path.basename(csv_rbh_output_filename))  # add csv file
        zf.write(csv_strict_output_filename, os.path.basename(csv_strict_output_filename))  # add csv file
        zf.write(csv_ns_output_filename, os.path.basename(csv_ns_output_filename))  # add csv file
    return zip_file


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


def file_len(fname):
    """Return the file length in lines."""
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def targz_folder(archive_name, folder):
    """
    Returns True after
    :param archive_name:
    :param folder:
    :return:
    """
    with tarfile.open(archive_name, "w:gz") as tar:
        tar.add(folder, arcname=os.path.basename(folder))
    return True


def cleanup(path, storage_folder, run_id):
    """
    Performs tar and gzip on sets of files produced by the program.
    Then deletes the files and folders.
    :param path:  # the run_folder
    :param storage_folder:  # the main folder, in which the entire run_folder will be stored
    :param run_id:  # the folder containing the first blast results
    :return:
    """
    # compress all files in path:
    # fasta_path
    path_archive = os.path.join(storage_folder, "{}.all.tar.gz".format(run_id))
    if targz_folder(path_archive, path):  # compress run_folder
        shutil.rmtree(path)  # delete run folder
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


def write_sort_command_script(filename_to_sort, sorted_filename):
    """Writing a sort uniq script to edit the gene csv file."""
    script_path = "/tmp/sort_script.sh"  # default script location
    with open(script_path, 'w') as script:
        script.write("#! /bin/tcsh\n")
        script.write("# The script is designed to run sort, uniq command from RecBlast\n")
        command_line = "sort {0} | uniq > {1}".format(filename_to_sort, sorted_filename)
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


def subset_db(tax_id, gi_file_path, db_path, big_db, run_anyway, DEBUG, debug, attempt_no=0):
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
                         target_db]  # TODO: test that blastdb_aliastool works for the user
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
re_sub = re.sub
upper = str.upper
lower = str.lower

