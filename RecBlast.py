#! /usr/bin/env python2

import argparse
# from RecBlastParams import *
from RecBlastUtils import *
import csv_transformer
import taxa_to_taxid
from uuid import uuid4
import part_one
import part_two
import part_three

# this will be the stand alone version of RecBlast for linux.

# flags:
# DEBUG = False
TAXA_ID_PROVIDED = False  # if the user provides a file with a list of taxa id.
GENE_CSV_PROVIDED = False
# If false, we will take a file of taxa names

# BLAST PARAMS
# defaults
E_VALUE_THRESH = 1e-7
IDENTITY_THRESHOLD = 37
COVERAGE_THRESHOLD = 50
MAX_TARGET_SEQS = '1000000'
BACK_E_VALUE_THRESH = E_VALUE_THRESH  # because it doesn't matter?
# BACK_E_VALUE_THRESH = 1e-15  # for the other direction
MAX_ATTEMPTS_TO_COMPLETE_REC_BLAST = 100
# CPU = 1

# fixed:
OUTFMT = '6 staxids sseqid pident qcovs evalue sscinames sblastnames'
ACCESSION_REGEX = re.compile(r'([A-Z0-9\._]+) ?')
ORIGINAL_ID = 1  # start part_two from 0. change this when you want to start from mid-file
APP_CONTACT_EMAIL = "recblast@gmail.com"
Entrez.email = APP_CONTACT_EMAIL

# comparison
TEXTUAL_MATCH = 0.4
TEXTUAL_SEQ_MATCH = 0.99

####################
#  parse arguments #
####################

parser = argparse.ArgumentParser()
parser.add_argument("origin_species", help="The species of origin (for which we start the blast)")
parser.add_argument("gene_file", help="A file containing a list of gene names to perform the reciprocal blast on.")
parser.add_argument("taxa_list_file", help="A file containing a list of taxa names to perform the reciprocal blast on. "
                                           "They must not include the original taxon!")
parser.add_argument("--output_path", help="A folder in which to keep all RecBlast output.")

parser.add_argument("--gene_csv", help="This flag means the gene file provided is already a CSV file containing the "
                                       "required genes as well as their description and uniprot id.",
                    action="store_true", default=GENE_CSV_PROVIDED)
parser.add_argument("--tax_ids", help="This flag means the taxa file provided is already a file containing the "
                                      "required taxa ids and not the taxa names.",
                    action="store_true", default=TAXA_ID_PROVIDED)

parser.add_argument("--run_name", help="The name the run will receive (will determine the folder names)")

parser.add_argument("--max_attempt_to_complete_recblast",
                    help="The maximum number of matches to perform the reciprocal blast on.",
                    default=MAX_ATTEMPTS_TO_COMPLETE_REC_BLAST)
# blast parameters
parser.add_argument("--num_threads", help="The number of threads (CPU) dedicated for parallel blast run.",
                    default=1, type=int)
parser.add_argument("--evalue", help="The e-value threshold for matches of the first blast.", default=E_VALUE_THRESH)
parser.add_argument("--evalue_back", help="The e-value threshold for matches of the second blast.",
                    default=BACK_E_VALUE_THRESH)
parser.add_argument("--identity", help="The minimum identity required for blast matches.",
                    default=IDENTITY_THRESHOLD)
parser.add_argument("--coverage", help="The minimum query and hit coverage required for blast matches.",
                    default=COVERAGE_THRESHOLD)
parser.add_argument("--max_seqs", help="The maximum number of sequences reported by blast.", default=MAX_TARGET_SEQS)
parser.add_argument("--db_first_run", help="The path to the BLASTP database for the first run.")
parser.add_argument("--target_db", help="The path to the BLASTP database for the second run.")

parser.add_argument("--string_similarity", help="The string similarity value for comparing the gene names/descriptions",
                    default=TEXTUAL_MATCH)

parser.add_argument("--run_even_if_no_db_found", help="Performs a heavy reciprocal blast. "
                                                      "Not recommended in most cases. See documentation for use cases.",
                    action="store_true")
parser.add_argument("--keep_files", help="Keeps intermediate files after completion", action="store_true",
                    default=False)
parser.add_argument("--debug", help="Adds debug prints in various stages of the run.", action="store_true",
                    default=False)

args = parser.parse_args()

# DEBUG flags
DEBUG = args.debug


def debug(s):
    return debug_s(s, DEBUG)


# making sure the files we received exist and are not empty:
if exists_not_empty(args.gene_file):
    debug("The gene file {} exists and not empty".format(args.gene_file))
else:
    print("gene_file {} does not exist or is empty!".format(args.gene_file))
    exit(1)

if exists_not_empty(args.taxa_list_file):
    debug("The taxa list file {} exists and not empty".format(args.taxa_list_file))
else:
    print("taxa_list_file {} does not exist or is empty!".format(args.taxa_list_file))
    exit(1)


# locating BLASTP path on your system
BLASTP_PATH = "Not valid"
try:
    BLASTP_PATH = subprocess.check_output(["which", "blastp"], universal_newlines=True).strip()
    debug("BLASTP found in {}".format(BLASTP_PATH))
except subprocess.CalledProcessError:
    print("No BLASTP found. Please check install blast properly or make sure it's in $PATH. Aborting.")
    exit(1)

CPU = args.num_threads

# parsing the genes_csv if it's a csv, and transforming it if it's a gene list file
if args.gene_csv:
    CSV_PATH = args.gene_file  # hope it's a good and valid file...
else:  # if the csv is not provided, create it from a gene file
    CSV_PATH = csv_transformer.gene_file_to_csv(args.gene_file)  # quits if no input could be converted

# script folder
SCRIPT_FOLDER = os.path.dirname(os.path.abspath(__file__))

# defining run folder
run_folder = os.getcwd()   # current folder
if args.run_name:
    run_name = args.run_name
    if args.output_path:
        run_folder = args.output_path  # assigning run_folder (default)
    else:
        run_folder = os.path.join(run_folder, run_name)
else:
    run_name = str(uuid4())  # randomly assigned run_name
    run_folder = os.path.join(run_folder, run_name)

create_folder_if_needed(run_folder)  # creating the folder

# creating the rest of the folders:
# folders:
FIRST_BLAST_FOLDER = os.path.join(run_folder, "first_blast")
create_folder_if_needed(FIRST_BLAST_FOLDER)  # creating the folder
SECOND_BLAST_FOLDER = os.path.join(run_folder, "second_blast")
create_folder_if_needed(SECOND_BLAST_FOLDER)  # creating the folder
FASTA_PATH = os.path.join(run_folder, "fasta_path")
create_folder_if_needed(FASTA_PATH)
FASTA_OUTPUT_FOLDER = os.path.join(run_folder, "fasta_output")
create_folder_if_needed(FASTA_OUTPUT_FOLDER)
CSV_OUTPUT_FILENAME = os.path.join(run_folder, "output_table.csv")
# Decide on taxa input:
TAX_DB = os.path.join(SCRIPT_FOLDER, "DB/taxdump/tax_names.txt")
# database location
DB_FOLDER = os.path.join(SCRIPT_FOLDER, "DB")


# parsing and creating taxa files and parameters:
tax_name_dict = taxa_to_taxid.create_tax_dict(TAX_DB)
tax_id_dict = dict((v, k) for k, v in tax_name_dict.iteritems())  # the reverse dict

# processing original species
if is_number(args.origin_species):  # already a tax_id
    ORG_TAX_ID = args.origin_species
    # write it to a file
    ORIGINAL_TAXA_FILE = os.path.join(run_folder, "original_tax_id.txt")
    with open(ORIGINAL_TAXA_FILE, 'w') as tax_file:
        tax_file.write(args.origin_species+"\n")
    ORIGIN_SPECIES = tax_id_dict[args.origin_species]
    debug("Tax id given: {0}, which is {1}, no need to convert. Saved original tax_id to {2}".format(
        args.origin_species, ORIGIN_SPECIES, ORIGINAL_TAXA_FILE))

else:  # it's a tax name
    ORIGIN_SPECIES = args.origin_species
    # convert it to tax id and write it to a file
    ORG_TAX_ID = tax_name_dict[ORIGIN_SPECIES]
    ORIGINAL_TAXA_FILE = os.path.join(run_folder, "original_tax_id.txt")
    with open(ORIGINAL_TAXA_FILE, 'w') as tax_file:
        tax_file.write(ORG_TAX_ID + "\n")
    debug("Tax name given: {0}, which is tax_id {1}, converted. Saved original tax_id to {2}".format(
        ORIGIN_SPECIES, ORG_TAX_ID, ORIGINAL_TAXA_FILE))

# working on the taxa list files
if args.tax_ids:  # already a tax_id file
    TAXA_LIST_FILE = args.taxa_list_file  # actually better to do them both
    # (TAXA_LIST_FILE, bad_tax_list) = taxa_to_taxid.convert_tax_to_taxid(tax_name_dict, args.taxa_list_file)
    debug("Taking fasta filters from pre made taxa list files. saved new list in: {}".format(TAXA_LIST_FILE))
else:  # only a taxa names file
    # converting taxa names list
    (TAXA_LIST_FILE, bad_tax_list) = taxa_to_taxid.convert_tax_to_taxid(tax_name_dict, args.taxa_list_file)
    if len(bad_tax_list) > 0:
        print("Bad taxa names found in the file provided:")
        print("\n".join(bad_tax_list))
        print("Ignoring them.")
    debug("Converted tax list to tax ID files, saved new list in: {}".format(TAXA_LIST_FILE))

# processing DB information:
#############################

# forward DB:
DB = "Not valid"  # default value. Doesn't matter at this point
if args.db_first_run:  # if the user provided a DB
    DB = args.db_first_run
else:  # find the nr DB yourself
    BLASTDB_PATH = ""
    try:
        BLASTDB_PATH = subprocess.check_output(["echo", "$BLASTDB"], universal_newlines=True).strip()
        if BLASTDB_PATH == "":  # depends on the system
            blastdb_exit()
    except subprocess.CalledProcessError:
        blastdb_exit()
    debug("Found $BLASTDB path in {}".format(BLASTDB_PATH))

    # check if nr exists
    if exists_not_empty(os.path.join(BLASTDB_PATH, "nr.00.phd")):
        DB = os.path.join(BLASTDB_PATH, "nr")
    else:
        print("Could not find the nr protein database in your $BLASTDB folder ({}).".format(BLASTDB_PATH))
        print("This means one of the following:\n"
              "1. You have the database but it's somewhere else.\n"
              "   If you do, please change the $BLASTDB variable to the nr location on your machine.\n"
              "2. You don't have the nr database at all.\n"
              "   If that's the case, please install it using this simple script: update_blastdb.pl\n"
              "   The script should already be installed as it comes with the Blast+ installation.\n"
              "Then run RecBlast again.")
        exit(1)

# target db (database of the original species)
if args.target_db:  # if the user provided
    TARGET_DB = args.target_db
else:
    TARGET_DB_FOLDER = os.path.join(DB_FOLDER, ORG_TAX_ID)  # this is where our db should be
    # check if it exists - if so use it as a db
    if os.path.exists(TARGET_DB_FOLDER):
        TARGET_DB = os.path.join(TARGET_DB_FOLDER, 'db')
        print "{} already has a local version of BLASTP DB!" .format(ORG_TAX_ID)
    else:  # if not create an alias
        print("No local version of {} database exists. Creating a subset now.".format(ORG_TAX_ID))
        gi_file = os.path.join(run_folder, "taxon_gi_file_list.txt")  # the path to the new file
        TARGET_DB = subset_db(ORG_TAX_ID, gi_file, DB_FOLDER, DB, args.run_even_if_no_db_found, DEBUG, debug)


#################
# Run main code #
#################

print "Welcome to RecBlast."
print "Run {0} started at: {1}".format(run_name, strftime('%H:%M:%S'))

# part 1:
print("starting to perform part_one.py")
id_dic, blast1_output_files = part_one.main(CSV_PATH, APP_CONTACT_EMAIL, run_folder, FASTA_PATH, FIRST_BLAST_FOLDER,
                                            FASTA_OUTPUT_FOLDER, BLASTP_PATH, DB, TAXA_LIST_FILE, OUTFMT,
                                            MAX_TARGET_SEQS, E_VALUE_THRESH, COVERAGE_THRESHOLD, CPU, DEBUG, debug)
print("BLASTP part 1 done!")  # should find a way to make sure the data is okay...
print("*******************")

# part 2:
second_blast_for_ids_dict, blast2_output_files, blast2_gene_id_paths = part_two.main(FIRST_BLAST_FOLDER,
                                                                                     SECOND_BLAST_FOLDER, ORIGINAL_ID,
                                                                                     E_VALUE_THRESH, IDENTITY_THRESHOLD,
                                                                                     COVERAGE_THRESHOLD,
                                                                                     ACCESSION_REGEX, run_folder,
                                                                                     BLASTP_PATH, TARGET_DB, OUTFMT,
                                                                                     MAX_TARGET_SEQS,
                                                                                     BACK_E_VALUE_THRESH, CPU,
                                                                                     ORIGINAL_TAXA_FILE, DEBUG, debug,
                                                                                     input_list=blast1_output_files)
print("BLASTP part 2 done!")
print("*******************")

# part 3:
if part_three.main(SECOND_BLAST_FOLDER, BACK_E_VALUE_THRESH, IDENTITY_THRESHOLD, COVERAGE_THRESHOLD, TEXTUAL_MATCH,
                   TEXTUAL_SEQ_MATCH, ORIGIN_SPECIES, ACCESSION_REGEX, run_folder, MAX_ATTEMPTS_TO_COMPLETE_REC_BLAST,
                   CSV_OUTPUT_FILENAME, FASTA_OUTPUT_FOLDER, DEBUG, debug, id_dic, second_blast_for_ids_dict,
                   blast2_gene_id_paths):
    print("part 3 done!")
    print("*******************")

# cleaning:
if not DEBUG and not args.keep_files:
    if cleanup(run_folder, FASTA_PATH, FIRST_BLAST_FOLDER, SECOND_BLAST_FOLDER):
        print("Files archived, compressed and cleaned.")

print("Program done.")


# TODO: test: csv conversion and tax_list conversion
# TODO create online parsing program
# TODO: create an online version
