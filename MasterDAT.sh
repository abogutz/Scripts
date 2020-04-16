#! /bin/bash
#$ -cwd
#$ -pe ncpus 4
#$ -l h_vmem=20G
#$ -m e
#$ -M tiffyyleung@gmail.com



###From tab-delimited file with SRACODE & NAME of dataset
###All the data are first downloaded, then trimmed (optional), and aligned 

#TODO: Provide full path to where Github Scripts directory is located 
#Ensure functions, MasterDAT.sh and config files are within same directory
SCRIPTS_DIR=/brcwork/lorincz_lab/tleung/Scripts

#############################################

#TODO: Ensure it's the correct config file for the server
source $SCRIPTS_DIR/Graham.config
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




