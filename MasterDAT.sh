#! /bin/bash


#############################################

#Choose the correct config file specific to the server you are currently using (see SCRIPTS_DIR for config files)
PROJECT_DIR="/project/def-mlorincz/"
source $PROJECT_DIR/scripts/MasterDAT/"ComputeCanada.config"
source $PROJECT_DIR/scripts/MasterDAT/"SRAtoBW_functions.sh"
# SHELL_SCRIPT=$SCRIPTS_DIR/$(basename $0) #TODO still needed?
export NCBI_API_KEY=e0a9648272551712054923626410c9135506 #This API key belongs to Aaron, and allows more unfettered access to NCBI servers

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




