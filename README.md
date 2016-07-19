# RecBlast - stand-alone version
RecBlast: Reciprocal Blast for the masses.

Welcome to RecBlast stand-alone.
RecBlast stand-alone is a local, downloadable version of RecBlast. RecBlast is a powerful tool to which performs large-scale orthology detection across multiple genes and multiple taxa (MGMT) by implementing fast reciprocal protein Blast.
RecBlast is available as an easy-to-use cloud-based web server, and can be accessed from [here](http://reciprocalblast.com).

The local version gives the user more freedom and flexibility and poses no limitations on the number of genes or taxa tested. Is is also completely open-source, and contributions are welcome.

### How does it work?
RecBlast performs reciprocal BLAST on several genes across several taxa. 
It performs BLASTP twice and filters the results according to quality parameters, providing a fast, easy way to discover orthologues across taxa.

A more detailed explanation of the pipeline, the analysis and the output can be found on the [detailed documentation page on the RecBlast website](http://reciprocalblast.com/explain).

### Hardware Requirements
RecBlast can run on any linux machine and doesn't have strict hardware requirements.
However, since it utilizes the BLAST+ multi-processing feature, there is a natural advantage to multi-core machines and clusters.
There's an additional advantage to computers with large RAM, especially when using multiple cores (the memory footprint is higher). 

### Download
In order to download the program, just clone the source from github using:
`git clone https://github.com/neuhofmo/RecBlast.git`

### Install
Since RecBlast is written in python, there's not much installation to do. Just cd into the RecBlast folder and run RecBlast.py.

### Requirements
    RecBlast can run on any linux machine.
    However, there are a few mandatory requirements for RecBlast to run properly:
#####1. BLASTP
RecBlast runs blastp from the BLAST+ suite, which means you need to have BLAST+ installed on your machine, and it should be from version 2.3.0 and above. 
If you don't have BLAST+ installed, you can get it from the NCBI FTP [here](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).
You can either download a binary copy or compile it from the source.
For any issue with installing BLAST+, you can go to their [user manual](https://www.ncbi.nlm.nih.gov/books/NBK279690/).
#####2. BLASTP has to be accessible through the "blastp" command.
If you install BLAST+ through apt-get (on Ubuntu) or compiled it from source, it most likely isn't a problem.
However, if you directly downloaded the binaries from the NCBI FTP, it might not be in your $PATH variable. 
In such case, you have several ways to solve this easy problem:
a. Copy blastp into one of your path folders (if you have `sudo` permissions):
`sudo cp <blastp> /usr/local/bin/blastp`
b. Add the current blastp folder to your path:
If you're using bash, just add the following line to ~/.bashrc:
`export PATH=$PATH:<your blastp folder>`
If you're using tcsh, just add the following line to ~/.tcshrc:
`setenv PATH $PATH:<your blastp folder>`
#####3. Local copy of the NR and taxdb databases
If you have BLAST+ already installed, please make sure you already have the NCBI nr database on your local machine, and that it is accessible by calling the $BLASTDB variable.
$BLASTDB should point to the folder containing those two databases. If it doesn't, make sure it does by setting the variable on startup:
`export BLASTDB=<your blast db folder>` on bash, and `setenv BLASTDB=<your blast db folder>` on tcsh (as described above).
If you don't have the NR database you can download it using the update_blastdb.pl command (comes with BLAST+).
`update_blastdb.pl --passive --decompress nr`
`update_blastdb.pl --passive --decompress taxdb`
You can also download those databases manually from the [NCBI FTP db folder](ftp://ftp.ncbi.nlm.nih.gov/blast/db/) and decompress them yourself. 
Don't forget to download the databases to the $BLASTDB folder!

#####4. Python 2.X (2.7 and above but not 3)
RecBlast currently runs on python 2 only.

#####5. Python packages
RecBlast depends on a few packages mentioned in requirements.txt
You can install them all using `pip`:
`sudo pip install -r requirements.txt`
and in case python2 is not your default:
`python2 -m pip install -r requirements.txt` (or any other pip command that works for you)
Or install them manually.
The packages needed are currently: [biopython](https://github.com/biopython/biopython.github.io/), uuid, argparse, requests, [mygene](https://pypi.python.org/pypi/mygene)
If you have any problem installing them, try following [this help page](http://stackoverflow.com/questions/26053982/error-setup-script-exited-with-error-command-x86-64-linux-gnu-gcc-failed-wit)

### Usage
RecBlast will create the output folder and save the output in your current directory.
In order to run RecBlast, just call RecBlast.py (`./RecBlast.py` or `python RecBlast.py` for example) with the following arguments:
`python RecBlast.py <origin_species> <gene_file> <taxa_list_file>`

Where:
origin_species:         The species of origin (for which we start the blast)
gene_file:              A file containing a list of gene names to perform the reciprocal blast on.
taxa_list_file:         A file containing a list of taxa names to perform the reciprocal blast on. They must not include the original taxon!

These should be enough to run RecBlast correctly. 
There is no real need for any other parameter to change, unless, for example, you have different threshold, you wish to use higher CPU or you already have your databases prepared in advance.

In addition, a list of optional arguments is described when using the help flag:
`python ./RecBlast.py --help`
Additional explanation of the effect of each argument can be found [here](http://reciprocalblast.com/documentation)

If you have any other questions, or any bugs to report, please contact us at recblast@gmail.com.

Good luck,

Efrat and Moran,
The RecBlast Team

