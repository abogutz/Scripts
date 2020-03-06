#! /bin/bash

###From tab-delimited file with SRACODE & NAME of dataset
###All the data are first downloaded, then trimmed (optional), and aligned 

source /media/barb/Heisenberg/Tiffany/Scripts/SRAtoBW_functions.sh
SHELL_SCRIPT=$(pwd)/$0

############### PIPELINE-SPECIFIC VARIABLES ###############



###############

parseOptions $@ #parse options from command line to this function
parallelRun
masterDownload
trimReads
masterAlign
collapseReplicates
masterTrackHub
