STEPS TO SET UP PIPELINE ON SERVER (To be performed by admin):

Download all required software for the pipeline, put in /projects/def-GROUPNAME/scripts/utilities
The ones listed are the ones required for Graham/Cedar.
entrez-direct: https://www.ncbi.nlm.nih.gov/books/NBK179288/
deeptools: https://deeptools.readthedocs.io/en/develop/content/installation.html
hicup (?)

Open the config file and complete all the USER ACTION.

Open the MasterDAT.sh and complete all the USER ACTION.


Also, edit the top of the script with the corresponding submission requirement. 
Fill in any <> spaces with your own information.

----> COMPUTE CANADA (CEDAR/GRAHAM)
#! /bin/bash
#SBATCH --account=<Group Name>            # required
#SBATCH --ntasks=8                        # number of MPI processes
#SBATCH --mem-per-cpu=12G                 # memory; default unit is megabytes
#SBATCH --time=01-12:00                   # time (DD-HH:MM)
#SBATCH --mail-user=<email address>
#SBATCH --mail-type=ALL

----> BRC
#! /bin/bash
#$ -cwd				# execute script from current directory
#$ -pe ncpus 4			# parallel environment
#$ -l h_vmem=15G		# memory
#$ -m e
#$ -M <email address>

STEPS TO RUN THE MASTERDAT PIPELINE
1. Create input file with datasets to be download separated by TAB.
Column 1 = SRACODE
Column 2 = NAME
Column 3 = STRAIN 1 (optional - for allele specific run)
Column 4 = STRAIN 2 (optional - for allele specific run)
EXAMPLE:
SRA1.1	NAME1_rep1
SRA1.2	NAME1_rep2
SRA2.1	NAME2_Rep1
SRA2.2	NAME2_Rep2

2. Run MasterDAT.sh in the current directory with options.
EXAMPLE:
../Scripts/MasterDAT.sh -i INPUT_FILE.txt -k

The shell script will separate the input file into different subsets and submit the different subsets onto the server itself.
The scripts waits 30 seconds between each submission.
