#! /bin/bash

###From tab-delimited file with SRACODE & NAME of dataset
###All the data are first downloaded, then trimmed (optional), and aligned 

source /media/barb/Heisenberg/Tiffany/Scripts/SRAtoBW_functions.sh

pushd $(dirname $0) > /dev/null
SHELL_SCRIPT=$(pwd -P)/$(basename $0)
popd > /dev/null

############### PIPELINE-SPECIFIC VARIABLES ###############



###############

parseOptions $@ #parse options from command line to this function
parallelRun
masterDownload
trimReads
masterAlign
collapseReplicates
masterTrackHub
