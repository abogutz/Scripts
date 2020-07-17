
SCRIPTS_DIR=/Users/jra/bin/scriptiff

source $SCRIPTS_DIR/jrambp.config
source $SCRIPTS_DIR/SRAtoBW_functions.sh
SHELL_SCRIPT=$SCRIPTS_DIR/$(basename $0)




parseOptions $@ # parse options from command line to this function
checkPseudogenome
parallelRun
setUp

masterDownload
masterAlign
postAlignmentCleanup

collapseSEPEbismark
collapseReplicatesBismark


masterTrackHub #convert bam to bigWigs and create TrackHub hierarchy




# OUTPUTS

# for each dataset:
	# PE and SE alignment statistics
	# PE and SE bams
	# PE and SE deduplicated bams
	# PE and SE CpG_report.txt
	# PE and SE coverage
	# PE and SE bisulphite conversion rate
	# SEPE mergeTwoStrands CpG report (for Keegan K's DMR finder)
	# SEPE CpGs covered (x1,2,3,4,5,10,15,20)
	# DNAme bigWigs
	
# for merged replicates:
	# SEPE, if exist, are all combined already
	# covered CpGs (x1,2,3,4,5,10,15,20)
	# CpG_report.txt
	# DNAme bigWigs