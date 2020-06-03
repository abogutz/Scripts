#! /bin/bash
#SBATCH --account=<Group Name>            # required
#SBATCH --ntasks=8                        # number of MPI processes
#SBATCH --mem-per-cpu=12G                 # memory; default unit is megabytes
#SBATCH --time=01-12:00                   # time (DD-HH:MM)
#SBATCH --mail-user=<email address>
#SBATCH --mail-type=ALL



###From tab-delimited file with SRACODE & NAME of dataset
###All the data are first downloaded, then trimmed (optional), and aligned 

### USER ACTION REQUIRED ###
#Provide full path to where Github Scripts directory is located
#Please ensure functions, MasterDAT.sh and config files are within same directory (don't move them!)
SCRIPTS_DIR=/project/def-mlorincz/tleung/Scripts

#############################################

### USER ACTION REQUIRED ###
#Choose the correct config file specific to the server you are currently using (see SCRIPTS_DIR for config files)
source $SCRIPTS_DIR/ComputeCanada.config
source $SCRIPTS_DIR/SRAtoBW_functions.sh
SHELL_SCRIPT=$SCRIPTS_DIR/$(basename $0)

############### PIPELINE ###############

parseOptions $@ #parse options from command line to this function
checkPseudogenome
parallelRun
setUp

masterDownload
trimReads
masterAlign

if $ALLELE_SPECIFIC; then
	ALLELE_RUN=true
	setPseudogenome #change reference genome of alignement to pseudogenome
	masterAlign #align fastq to pseudogenome
	unpackAllelic $HAPLO_1 #unpack reads from the aligned files into two different files to look at allele specific
	unpackAllelic $HAPLO_2
fi

collapseReplicates
projectAllelic
masterTrackHub
removeFASTQ




