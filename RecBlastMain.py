#! /usr/bin/env python2
# wrapper for the other scripts

from RecBlastParams import *
import part_one
import part_two
import part_three

print "Welcome to RecBlast."
# print "Run {0} (id: {1}) started at: {2}".format(run_name, run_id, strftime('%H:%M:%S'))

# part 1:
print("starting to perform part_one.py")
id_dic, blast1_output_files = part_one.main()  # get job_num_list and check jobs between each script

print("BLASTP part 1 done!")  # should find a way to make sure the data is okay...
print("*******************")

second_blast_for_ids_dict, blast2_output_files, blast2_gene_id_paths = part_two.main(input_list=blast1_output_files)
print("BLASTP part 2 done!")  # should find a way to make sure the data is okay...
print("*******************")


# part 3
if part_three.main(id_dic, second_blast_for_ids_dict, blast2_gene_id_paths):
    print("part 3 done!")
    print("*******************")

# cleaning:
if not DEBUG and cleanup(run_folder, FASTA_PATH, FIRST_BLAST_FOLDER, SECOND_BLAST_FOLDER):
        print("Files archived, compressed and cleaned.")


print("Program done.")


# added the qsub -q low.q
# replace gis with accession.version
# debug 3: only 1 organism
# add unique hash to every run
# get taxa_filter to gzip the full files...
# Change update_match_results format
# add email support
# add email address for our ncbi connection
# test email module...
# change pickle folder
# edit filter taxa to work with a list
# integrate additional scripts to run before the gene_list_to_csv program
# add grep -w -f ../../../taxa_human.txt secondblast_fasta_for_id_2_index_500_full.txt | wc -l
# add filter taxa to the run itself...
# make sure we have the taxa_id file!!!
# add fasta update_match_results  - waiting for Efrat
# debug part 2 - the counters keep running after moving to the next id - hopefully solved by changing to global
# added rbh update_match_results to regular csv update_match_results (untested!)
# debug 3: [] incrementation, no success.
# run name in addition to run_id
# delete intermediate files
# clean code
# accept user parameters somehow
# change defaults - depending on the machine it's running on
# test runs
# try part 2,3
# edit parameters
# add db logic
# add usage and argparse (in the other main)

