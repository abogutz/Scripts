#! /bin/bash
#SBATCH --account=<Group Name>            # required
#SBATCH --ntasks=8                        # number of MPI processes
#SBATCH --mem-per-cpu=12G                 # memory; default unit is megabytes
#SBATCH --time=01-12:00                   # time (DD-HH:MM)
#SBATCH --mail-user=<email address>
#SBATCH --mail-type=ALL



# Modify these to whatever you want
INPUT_FILE=
GENOME=

# Many other options exist for MasterDAT.sh, run with -h option to see all of them

/project/def-mlorincz/scripts/MasterDAT/MasterDAT.sh -i $INPUT_FILE -g $GENOME