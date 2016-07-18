#! /usr/bin/env python2

from RecBlastUtils import *

# flags:
DEBUG = True
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
CPU = 1

# fixed:
OUTFMT = '6 staxids sseqid pident qcovs evalue sscinames sblastnames'
ACCESSION_REGEX = re.compile(r'([A-Z0-9\._]+) ?')
ORIGINAL_ID = 1  # start part_two from 0. change this when you want to start from mid-file
APP_CONTACT_EMAIL = "recblast@gmail.com"
Entrez.email = APP_CONTACT_EMAIL

# comparison
TEXTUAL_MATCH = 0.4
TEXTUAL_SEQ_MATCH = 0.99

# GENES_FILE = "/export/groups/igv/moranne/RecBlast/815_check_list.csv"  # our big run
# GENES_FILE = "/export/groups/igv/moranne/RecBlast/21_check_list.csv"  # our big run  # short real list from 815 (...81ac)
# GENES_FILE = "/export/groups/igv/moranne/RecBlast/11_control_list.csv"  # short control list from 447 (...497f)

# TAXA_LIST_FILE = "full_taxa_list.txt"  # doing another one for testing
# ORIGINAL_TAXA_FILE = "taxa_human.txt"  #

