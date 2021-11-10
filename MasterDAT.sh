#! /bin/bash


#############################################

#Choose the correct config file specific to the server you are currently using (see SCRIPTS_DIR for config files)
source "ComputeCanada.config"
source "SRAtoBW_functions.sh"
# SHELL_SCRIPT=$SCRIPTS_DIR/$(basename $0) #TODO still needed?

############### PIPELINE ###############

parseOptions $@ #parse options from command line to this function
loadModules
checkPseudogenome
#parallelRun
setUp

masterDownload
trimReads
masterAlign

if $ALLELE_SPECIFIC; then
	ALLELE_RUN=true
	setPseudogenome #change reference genome of alignment to pseudogenome
	masterAlign #align fastq to pseudogenome
	unpackAllelic $HAPLO_1 #unpack reads from the aligned files into two different files to look at allele specific
	unpackAllelic $HAPLO_2
fi

collapseReplicates
projectAllelic
masterTrackHub
removeFASTQ




