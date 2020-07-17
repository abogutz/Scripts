#!/bin/bash

#SBATCH --account=def-mlorincz           # required
#SBATCH --nodes=1                        # cores distributed in # of nodes
#SBATCH --ntasks-per-node=8              # number of MPI processes
#SBATCH --mem=84G                        # total memory; default unit is megabytes
#SBATCH --time=00-1:00                  # time (DD-HH:MM)

# High-level script that calls functions from SRAtoBW_functions.sh
# Takes as input a tab-delimited file with SRACODE & NAME of datasets (run MasterDAT.sh -h for help)
# By default, data are first downloaded, trimmed (optional), aligned and converted to normalized bigWigs


### USER ACTION REQUIRED ###
#Provide full path to where Github Scripts directory is located
#Please ensure functions, MasterDAT.sh and config files are within same directory (don't move them!)
SCRIPTS_DIR=/Users/jra/bin/scriptiff

### USER ACTION REQUIRED ###
#Choose the correct config file specific to the server you are currently using (see SCRIPTS_DIR for config files)
source $SCRIPTS_DIR/jrambp.config
source $SCRIPTS_DIR/SRAtoBW_functions.sh
SHELL_SCRIPT=$SCRIPTS_DIR/$(basename $0)



# SET MAPQ HERE?


############### PIPELINE ###############

parseOptions $@ #parse options from command line to this function
checkPseudogenome
parallelRun
setUp

masterDownload
trimReads
masterAlign
collapseReplicates #combine bam files generated from technical or biological duplicates (end in _rep1 _rep2 or _Rep1 _Rep2)

masterTrackHub #convert bam to bigWigs and create TrackHub hierarchy
# removeFASTQ




