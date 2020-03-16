  #! /bin/bash
#$ -cwd
#$ -pe ncpus 4
#$ -l h_vmem=10G
#$ -m e
#$ -M tiffyyleung@gmail.com



###From tab-delimited file with SRACODE & NAME of dataset
###All the data are first downloaded, then trimmed (optional), and aligned 

#############################################
########## BEFORE USING THE SCRIPT ##########
#############################################

#TODO: provide full path to directory holding functions & config files
#Ensure function and config files are within same directory (sourcing will not work otherwise)
#FUNCTIONS_DIR=/media/barb/Heisenberg/Tiffany/Scripts/
FUNCTIONS_DIR=/brcwork/lorincz_lab/tleung/Scripts/


#TODO: alter any system specific variables and tools path through config file

#############################################


source $FUNCTIONS_DIR/SRAtoBW_functions.sh
#source $FUNCTIONS_DIR/MEA_functions.sh

# The following 3 lines of code is to obtain the full path of this script for parallel submission 
pushd $(dirname $0) > /dev/null
SHELL_SCRIPT=$(pwd -P)/$(basename $0)
popd > /dev/null
 


############### PIPELINE-SPECIFIC VARIABLES ###############



############### PIPELINE FUNCTIONS ###############

parseOptions $@ #parse options from command line to this function
parallelRun
masterDownload
trimReads
masterAlign
collapseReplicates

if $ALLELE_SPECIFIC; then
setPseudogenome #change reference genome of alignement to pseudogenome
masterAlign #align fastq to pseudogenome
collapseReplicates
allespecUnpack #unpack reads from the aligned files into two different files to look at allele specific 
fi

masterTrackHub
cleanFASTQ



