  #! /bin/bash
#$ -cwd
#$ -pe ncpus 3
#$ -l h_vmem=12G
#$ -m e
#$ -M tiffyyleung@gmail.com



###From tab-delimited file with SRACODE & NAME of dataset
###All the data are first downloaded, then trimmed (optional), and aligned 

#############################################
########## BEFORE USING THE SCRIPT ##########
#############################################

#TODO: alter any system specific variables and tools path through config file
#TODO:Ensure functions, MasterDAT.sh and config files are within same directory (sourcing will not work otherwise)
if [[ -z $FUNCTIONS_DIR ]]; then #if want to use functions by themselves
  pushd $(dirname $0) > /dev/null
  FUNCTIONS_DIR=$(pwd -P)
  popd > /dev/null
fi

#############################################

SHELL_SCRIPT=$FUNCTIONS_DIR/$(basename $0)
source $FUNCTIONS_DIR/SRAtoBW_functions.sh

############### PIPELINE ###############

parseOptions $@ #parse options from command line to this function
checkPseudogenome
parallelRun
checkDependencies

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




