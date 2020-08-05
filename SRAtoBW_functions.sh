#! /bin/bash

#List of functions used in various pipelines for download, alignments and creating track hubs

# Pipelines that use these functions:
# MasterDAT.sh
# ...
# ...

#Troubleshooting
#Only change if you ever plan on using these functions directly in the command line (otherwise ignore)
#NO user action required here
#Ensure correct config file if using
if [[ -z $SCRIPTS_DIR ]]; then
	pushd $(dirname $0) > /dev/null
	SCRIPTS_DIR=$(pwd -P)
	popd > /dev/null
	source $SCRIPTS_DIR/jrampb.config
fi

# Non-system specific variables
CURRENT_DIRECTORY=$(pwd)


ALLELE_SPECIFIC=false
ALLELE_RUN=false
BAM_INPUT=false
BAM_REALIGNMENT=false
FASTQ_INPUT=false
FASTQ_ONLY=false
KEEP_FASTQ=true
KEEP_REPLICATES=true
LOCAL_RUN=true
SEP_PARA=false
STRANDED_ALLELIC=false
TRIM_READ=true		# trim by default
USE_BOWTIE=true 	# default is bowtie2

CODE_ARRAY=""
BIN_SIZE=1
FLAG=1540 #read unmapped, read fails platform/vendor quality checks, read is PCR or optical duplicate
GENOME_BUILD="mm10"
MIN_MAPQ=1
NORMALIZE="CPM"
SMOOTH_WIN=0

DEPENDENCIES=($ESEARCH $EFETCH $FASTERQDUMP "$TRIMMOMATIC" $STAR $BISMARK $BOWTIE2 $SAMTOOLS "$PICARD" awk $BAM2FASTQ $BEDGRAPHTOBW $BAMCOVERAGE)

# Help Menu
OPTIONS="hi:ab:B:d:Df:Fg:kLm:M:n:N:ors:t:Tux:X"

HELP="USAGE:\t $(basename $0) [OPTIONS] -h for help"

HELP_FULL="\n$HELP\n
\nThis script has been designed to streamline and standardize the procurement of data from the SRA database, as well as alignment and visualization of all data. Given an input file listing the SRA files to be downloaded, it will automatically download them, decompress, rename according to the users desires, align to the genome of choice, collapse replicates and finally generate a UCSC-compatible track hub including bigwig files. Users can opt to use some or all of these features depending on their input or output needs.\n
\nIt is generally best to run this script in a clean folder, as it performs many cleanup steps on both the current folder and ALL subdirectories. For instance, any folders that include \"STAR\" in the name will be removed, as the STAR program leaves temporary aligns behind when finished.\n
\nAlignment:\tChIPseq - bwa mem\n\t\tRNAseq - STAR\n\t\tBSSeq/PBAT/RRBS - bismark\n

\nPlease include the following information in your filenames:\n<RNAseq/ChIPseq/BSseq>_<study>_<replicate>
\nExample input.txt file:\nGSM3597247\tC57BL6J_oocyte_RNApolII_ChIPseq_Bogutz2019_rep1\nExample run command:\n~/bin/Scripts/MasterDAT.sh -i input.txt\n\n

\n\nOPTIONS:
-h\tPrints help page.\n\t
-i\tSRA Input File. Must be in the format Accession(tab)DesiredName.\n\t\tDesiredName must be formatted specifically. A final identifier \n\t\tmust be present at the end of the name by which data should \n\t\tbe grouped (ex. Foo2018, etc.). The aligned files will be in a \n\t\tdirectory of this name. If there are replicates, the name will\n\t\tbe followed by an underscore and \"RepX\" - this will be removed \n\t\tfrom the final name upon collapsing replicates. Data type \n\t\tshould be included somewhere in the name; if not, it will be\n\t\tprepended to the name as 'RNAseq' or 'BSSeq' or 'ChIPseq'.\n\t
-a\tAllele-specific alignment using MEA.\n\t
-b\tBAM Inputs. Will use any already aligned .bam files in\n\t\tsubdirectories to generate trackhub.\n\t
-B\tBAM Input with realignment. Will extract reads from .bam files\n\t\tin listed directory and realign to genome of choice.\n\t
-d\tTemporary directory. Useful for solid-state drives etc.\n\t\tPlease change the default using BRC.config.\n\t
-D\tCheck Dependencies and exit.\n\t
-f\tInput .fastq files. Provide folder name in which fastq files\n\t\tare located (files must end in .fastq.gz).\n\t\tData type must be included in name. Accepted: RNA, ChIP, BSSeq,\n\t\tPBAT, RRBS, ATAC, DNAse, HiC\n\t
-F\tOnly output .fastq files.\n\t
-g\tGenome build for alignment. Allowable: mm9, mm10, rn5, rn6,\n\t\toryCun2, mesAur1, or hg19. Default=mm10\n\t
-k\tKeep .fastq files when done.\n\t
-L\tRunning script on a local computer rather than a server.\n\t
-m\tMemory to give each thread (Format=XG/M).\n\t
-M\tMinimum mapping quality for bigwig generation. Default=5\n\t
-n\tBin size for bigwig generation. Larger bins to smooth noisy\n\t\tdata. Default=1\n\t
-N\tNormalization method for bigwigs. Accepted: CPM, RPKM\n\t\t(Default=CPM)\n\t
-r\tKeep replicates after collapsing. Default=false.\n\t
-s\tObtained stranded RPM tracks for allele-specific runs.\n\t\tDefault=false\n\t
-t\tNumber of Threads to use. Default=6 (Check in config file)\n\t
-T\tTrim .fastq files using Trimmomatic.\n\t
-u\tChanging ChIPseq aligner to bowtie2. Default=BWA\n\t
-x\tPrefix of fastq.gz files (used with -f)\n\t
-w\tSmoothing window. Will smooth bigwigs in a rolling window of\n\t\tthis size. Default=0\n\t"



########################################
#	      FUNCTIONS		       #
########################################

### Change/parse options/variables based on user-input parameters passed from command-line
function parseOptions () {
	if ( ! getopts $OPTIONS opt); then
		echo -e $HELP
		exit 1
	fi

	while getopts $OPTIONS opt; do
		case $opt in
			h) #open the HELP menu
				echo -e $HELP_FULL | fold -s
				exit
				;;
			i) #set input file
				SEP_PARA=true
				PASS_ARG=$@
				FILENAME=${OPTARG}
				;;
			a)
				ALLELE_SPECIFIC=true
				;;
			b) #bam input for trackhub
				BAM_INPUT=true
				if $FASTQ_ONLY ; then
					echo -e "ERROR:\tIncompatible options. Cannot generate .fastq files from .bam files."
					exit 1
				fi	
				;;
			B) #bam input for realignmnet
				BAM_REALIGNMENT=true
				BAM_REALIGNMENT_DIRECTORY=${OPTARG}
				;;
			d)
				TEMP_DIR=${OPTARG}
				;;
			D)
				DEPEND=1
				;;
			f) #use existing fastq files for downstream modules in pipeline (will skip download)
				FASTQ_INPUT=true
				if $FASTQ_ONLY ; then
					echo -e "ERROR:\tFiles already in Fastq format."
					exit 1
				fi
				FASTQ_DIRECTORY=${OPTARG}
				;;
			F) #stop pipeline after downloading
				FASTQ_ONLY=true
				if $BAM_INPUT ; then
					echo -e "ERROR:\tIncompatible options. Cannot generate .fastq files from .bam files."
					exit 1
				fi
				;;
			g)
				GENOME_BUILD=${OPTARG}
				;;
			k) #keep fastq files after alignment
				KEEP_FASTQ=true
				;;
			L)
				LOCAL_RUN=true
				;;
			m)
				THREAD_MEM=${OPTARG}
				;;
			M)
				MIN_MAPQ=${OPTARG}
				;;
			n)
				BIN_SIZE=${OPTARG}
				;;
			N)
				NORMALIZE=${OPTARG}
				if [[ $NORMALIZE != " " ]] && [[ $NORMALIZE != "RPKM" ]] ; then
					echo "Not an acceptable normalization method."
					exit 1
				fi
				;;
			r)
				KEEP_REPLICATES=true
				;;
			s)
				STRANDED_ALLELIC=true
				;;
			t)
				RUN_THREAD=${OPTARG}
				;;
			T) #trim fastq using TRIMMOMATIC and default parameters
				TRIM_READ=true
				;;
			u)
				USE_BOWTIE=true
				;;
			w)
				SMOOTH_WIN=${OPTARG}
				;;
			x) #used in conjuction with -f. Run pipeline only on files with prefix SEARCH_KEY
			   #is equivalent to USER_INPUT_NAME when using -i input_text_file.txt, WITHOUT "_rep"
				SEP_PARA=false
				SEARCH_KEY=${OPTARG}
				;;
			X)
				CODE_ARRAY=$(echo $@ | sed 's/.*-X //')
				;;
			\?)
				echo -e "\n###############\nERROR: Invalid Option! \nTry '$(basename $0) -h' for help.\n###############" >&2
				exit 1
				;;
		esac
	done

	setGenome $GENOME_BUILD
}

### Set up log file for parallel runs
function setUpLogFile () {
	if [[ $SEP_PARA == false ]]; then
		LOG_FILE=$CURRENT_DIRECTORY/$SEARCH_KEY"_"$(date '+%y-%m-%d')"_log.txt"
		checkDependencies
		printProgress "[setUpLogFile] Script started at [$(date)]"
		printProgress "[setUpLogFile] Search key for the set: $SEARCH_KEY"
		printProgress "[setUpLogFile] SRA array: $CODE_ARRAY"
		printProgress "[setGenome] Genome used for data is $GENOME_BUILD"
		mkdir -p $OUTPUT_BAM_FOLDER
		mkdir -p $OUTPUT_BIGWIG_FOLDER
		mkdir -p $OUTPUT_STATS_FOLDER
		mkdir -p $OUTPUT_CPG_REPORT_FOLDER
	fi
}

function printProgress () {
	echo -e $1 | tee -a $LOG_FILE #using tee will also show in stdout
}

function checkFileExists () {
	if [[ ! -f $1 ]]; then
		echo "ERROR: File $1 does not exist"
		exit 1
	fi
}

function checkBamFileExists () {
	if [[ ! -f $1 || $($SAMTOOLS view $1 | wc -l | cut -d ' ' -f4) == 0 ]]; then
		echo "ERROR: Bam file $1 does not exist or contains no aligned reads"
		exit 1
	fi
}

### Create neccessary files for reference genome
function setGenome () {
	if [[ $FASTQ_ONLY == false ]] ; then
		mkdir -p $TRACK_HUB_DIR
		printf "hub <HubNameWithoutSpace>\nshortLabel <max 17 char, display on side>\nlongLabel Hub to display <fill> data at UCSC\ngenomesFile genomes.txt\nemail <email-optional>" > ./$TRACK_HUB_DIR/hub.txt
	fi

### USER ACTION REQUIRED ###
# Please specify the exact location of your existing reference genomes and their indices
# Otherwise, ignore (right?)
	MOUSE=""
	RAT="."

	case $1 in
		"mm10")
			GENOME_DIR=$GENOME_DIR/$MOUSE/$1/
			if [[ $FASTQ_ONLY == false ]] ; then
#				printf "chr5\t142904365\t142904366" > Actb.bed
				TRACK_FOLDER=$TRACK_HUB_DIR/$1
				mkdir -p $TRACK_FOLDER
				TRACKDB=$TRACK_FOLDER/"trackDb.txt"
				printf "genome mm10\ntrackDb mm10/trackDb.txt" > $TRACK_HUB_DIR/genomes.txt
			fi
			;;
		"rn6")
			GENOME_DIR=$GENOME_DIR/$RAT/$1/
			if [[ $FASTQ_ONLY == false ]] ; then
				printf "chr12\t13718023\t13718024" > Actb.bed
				TRACK_FOLDER=$TRACK_HUB_DIR/$1
				mkdir -p $TRACK_FOLDER
				TRACKDB=$TRACK_FOLDER/"trackDb.txt"
				printf "genome rn6\ntrackDb rn6/trackDb.txt" > $TRACK_HUB_DIR/genomes.txt
			fi
			;;
		*)
			echo -e "ERROR: \t$1 is not a valid genome build. Enter -h for help."
			exit 1
			;;
	esac

	CHROM_SIZES=$GENOME_DIR/$1".sizes"
	GENOME_FILE=$GENOME_DIR/$1".fa"
	STAR_GENOME_DIR=$GENOME_DIR/$1"-STAR"
	BISMARK_GENOME_DIR=$GENOME_DIR
	BOWTIE2_INDEXES=$GENOME_DIR/$GENOME_BUILD

	DIPLOID_GENOME_DIR=$GENOME_DIR/diploid
	HAPLOID_GENOME_DIR=$GENOME_DIR/haploid

	checkFileExists $CHROM_SIZES
	checkFileExists $GENOME_FILE
}


### Create an array of SRACODE passed listed in the tab-delimited file
function createCodeArray () {
	while read n; do #read command will read a line of and split them into words (2 words in files)
		declare -a m=($n) #-a for an array m to be stored based on the line (SRACODE, NAME)
		CODE_ARRAY=$CODE_ARRAY" "${m[0]} #add only SRACODE to this list
	done < $INPUT_SAMPLE_LIST_FILE #feeding FILE as stdin to this while read
}


### Create subset of SRACODE array to be used in parallel running
function parallelRun () {
	
	# SEP_PARA is used to determine if the run is already in the parallel run or needs to be separated.
	# -i (include input file), SEP_PARA will be true. submit a new MasterDat with -X option to force it to be FALSE
	# i.e. the first MasterDat simply calls other MasterDat scripts. 
	if $SEP_PARA; then
		echo -e "Separating input file into subsets for parallel runs..."

		INPUT_SAMPLE_LIST_FILE=$CURRENT_DIRECTORY/$FILENAME
		createCodeArray
		CURRENT_SET=""

		for SRACODE in $CODE_ARRAY; do
			USER_INPUT_NAME=$(grep -e $SRACODE $INPUT_SAMPLE_LIST_FILE | cut -f2)

			if [[ $USER_INPUT_NAME == *_[Rr]ep* ]] ; then
				if [[ $CURRENT_SET != ${USER_INPUT_NAME//_[Rr]ep*/}* ]]; then
					CURRENT_SET=${USER_INPUT_NAME//_[Rr]ep*/}
					declare -a SUB_ARRAY=$(grep -e $CURRENT_SET $INPUT_SAMPLE_LIST_FILE | cut -f1)
					echo "calling "$(basename $SHELL_SCRIPT) "on" $CURRENT_SET

					if $LOCAL_RUN; then
						$SHELL_SCRIPT $PASS_ARG -x $CURRENT_SET -X $SUB_ARRAY &
						wait $! #wait for the script above to finish running before moving onto the next set (avoid overload)
					else
						echo "Submitting on server"
						DATE=$(date '+%y-%m-%d')
						$SERVER_SUBMIT "MasterDAT_"$DATE"_"$CURRENT_SET $SHELL_SCRIPT $PASS_ARG -x $CURRENT_SET -X $SUB_ARRAY
						sleep 30 #pause for 30 secs before running next code b/c fetching data takes some time
					fi
				fi

			else
				declare -a SUB_ARRAY=$SRACODE
				echo "calling "$(basename $SHELL_SCRIPT) "on" $USER_INPUT_NAME

				if $LOCAL_RUN ; then
					$SHELL_SCRIPT $PASS_ARG -x $USER_INPUT_NAME -X $SUB_ARRAY
					wait $!
				else
					echo "Submitting on server"
					DATE=$(date '+%y-%m-%d')
					$SERVER_SUBMIT "MasterDAT_"$DATE"_"$USER_INPUT_NAME $SHELL_SCRIPT $PASS_ARG -x $USER_INPUT_NAME -X $SUB_ARRAY
					sleep 30
				fi
			fi
			
		done
		exit #exit the script b/c don't want to run the rest of the code on every single SRACODE again
	fi
}


function setupVariables () {
	
	# These variables are defined in the function parallelRun
	# USER_INPUT_NAME=$(grep -e $SRACODE $INPUT_SAMPLE_LIST_FILE | cut -f2)
	# SEARCH_KEY=${USER_INPUT_NAME//_[Rr]ep*/}
	
	FULL_FILENAME=$1
	if [[ $FULL_FILENAME == *"_2.fastq.gz" ]] ; then
		continue 												#dont process 2nd read, continue to next iteration
	fi

	if [[ $FULL_FILENAME == *"_1.fastq.gz" ]] ; then
		PAIRED_END=true
		FASTQ_PATH=${FULL_FILENAME//_1.fastq.gz/} 				#removes everything that's after // - leaving path to directory
		FILE_FASTQ1=$FASTQ_PATH"_1.fastq.gz"
		FILE_FASTQ2=$FASTQ_PATH"_2.fastq.gz"

	else # SE sequencing
		PAIRED_END=false
		FASTQ_PATH=${FULL_FILENAME//.fastq.gz}
		FILE_FASTQ=$FASTQ_PATH".fastq.gz"
		USER_INPUT_NAME=$(basename $FASTQ_PATH)
	fi
	
	# must set this for runs that do not input an SRA+filename file
	USER_INPUT_NAME=$(basename $FASTQ_PATH)
	FILE_RAW_BAM=$USER_INPUT_NAME"_raw.bam"
	FILE_BAM=$USER_INPUT_NAME".bam"
	# This will become the new working directory
	
	
	# reset TMP_DIR variable to avoid compounding effect when FASTQ_INPUT=true... this is so bad. but it works..
	# maybe call it WORKING_DIR? and replcae all instances of TEMP_DIR with it?
	TEMP_DIR="/Users/jra/bioinf/scratch"
	NEW_TEMP_DIR_NAME=${SEARCH_KEY##*_} 						# removes everything before the last _ (leaving grouping identifier)
	TEMP_DIR=$TEMP_DIR/$NEW_TEMP_DIR_NAME
	mkdir -p $TEMP_DIR
	
	# Determine library type (pretty unnecessary... can just check file names directly!)
	if [[ $FULL_FILENAME == *"C[Hh]IP"* ]] ; then
		LIBRARY_TYPE="ChIPseq"
	elif [[ $FULL_FILENAME == *"RNA"* ]] ; then
		LIBRARY_TYPE="RNAseq"
	elif [[ $FULL_FILENAME == *"ATAC"* ]] ; then
		LIBRARY_TYPE="ChIPseq"
	elif [[ $FULL_FILENAME == *"DNA[Ss][Ee]"* ]] ; then
		LIBRARY_TYPE="ChIPseq"
	elif [[ $FULL_FILENAME == *"H[Ii]C"* ]] ; then
		LIBRARY_TYPE="HiC"
	elif [[ $FULL_FILENAME == *"RRBS"* ]] || \
		[[ $FULL_FILENAME == *"WGBS"* ]] || \
		[[ $FULL_FILENAME == *"DNA[Mm][Ee]"* ]] || \
		[[ $FULL_FILENAME == *"PBAT"* ]] || \
		[[ $FULL_FILENAME == *"[Mm]ethyl[Cc]"* ]] ; then
		LIBRARY_TYPE="Bisseq"
	else
		printProgress "[setupVariables] Unknown library type. Please provide a valid library type (see help menu) in the filename."
	fi
		
	echo -e "Testing: \nFull Filename: $FULL_FILENAME\nFASTQ_PATH: $FASTQ_PATH\nUSER_INPUT_NAME: $USER_INPUT_NAME \nSEARCH_KEY: $SEARCH_KEY \nPAIRED_END: $PAIRED_END \nFILE_BAM: $FILE_BAM\nTEMP_DIR: $TEMP_DIR\nLIBRARY_TYPE: $LIBRARY_TYPE"
}




### Downloading files specified from the tab-delimited file to fastq files
### When calling function by itself, $1=SRA input file
function masterDownload () {
	cd $TEMP_DIR
	mkdir -p $FASTQ_DIRECTORY
	
	if $BAM_INPUT || $FASTQ_INPUT; then #if input is either BAM/FASTQ then no download is needed - exit masterDownload function
		return 0
	fi
	
	if $BAM_REALIGNMENT; then
		obtainFastqFromBAM
		return 0 #after extracting fastq from BAM, can exit masterDownload function
	fi

	INPUT_SAMPLE_LIST_FILE=$CURRENT_DIRECTORY/$FILENAME
	
	printProgress "[masterDownload] Started at [$(date)]"
	for code in $CODE_ARRAY; do
		downloadReadsJRA $code
	done

	printProgress "[masterDownload] All fastq files for $SEARCH_KEY are downloaded [$(date)]"

	if $FASTQ_ONLY; then
		printProgress "[masterDownload] Fastq files only. Exit script."
		exit 0 #exit the whole script b/c only requires fastq files
	fi
}



function downloadReadsJRA () {
	SRACODE=$1
	USER_INPUT_NAME=$(grep -e $SRACODE $INPUT_SAMPLE_LIST_FILE | cut -f2)

	printProgress "[masterDownload one step JRA test] Downloading $SRACODE to temporary folder $TEMP_DIR... at [$(date)]"
	DL=$($ESEARCH -db sra -query $SRACODE | $EFETCH --format runinfo | cut -d ',' -f 1 | grep -v "Run")

	# DETERMINE PAIRED-ENDEDNESS
	COUNTING_DIR=$TEMP_DIR/$USER_INPUT_NAME"_counting_tmp"
	$FASTQDUMP $DL -X 1 -O $COUNTING_DIR --split-files
	if [[ ! -z `ls $COUNTING_DIR/*_2.fastq` ]]; then
		READ_1_COUNT=$(ls -l  $COUNTING_DIR/*_1.fastq | wc -l )
		READ_2_COUNT=$(ls -l  $COUNTING_DIR/*_2.fastq | wc -l )
		#for example GSM1845258 contains both SE And PE libraries. Complete horse shit.
		if [[ $READ_1_COUNT != $READ_2_COUNT ]]; then
			printProgress "Single-end runs detected alongside paired-end in $USER_INPUT_NAME. We recommend downloading runs individually. Passing to next entry"
			return 0 # leave behind temporary counting directory for troubleshooting
		else
			PAIRED_END=true
			printProgress "[masterDownload] Data are paired-end."
		fi
	else
		PAIRED_END=false
		printProgress "[masterDownload] Data are single-end."
	fi
	rm -r $COUNTING_DIR

	# DOWNLOAD
	# MAKE SURE DOWNLOADING IN PARALLEL IS NOT BUGGED USING THIS STRATEGY
	DOWNLOAD_DIR_TMP=$TEMP_DIR/"$USER_INPUT_NAME"_download_tmp
	if $PAIRED_END; then
		$FASTERQDUMP -O $DOWNLOAD_DIR_TMP -S -p $DL
	else
		$FASTERQDUMP -O $DOWNLOAD_DIR_TMP -p $DL
	fi

	# COMPRESS FASTQ FILES: THIS SHOULD BE A SEPARATE FUNCTION!
	printProgress "[masterDownload] Compressing fastq files..." #simultaneously zip all the dump fastq files for that read
	for DL_FASTQ in $DOWNLOAD_DIR_TMP/*fastq; do
		gzip $DL_FASTQ &
		COMPRESS_PID_ARRAY[$COMPRESS_COUNTER]=$!
		((COMPRESS_COUNTER++))
	done
	for gzipid in ${COMPRESS_PID_ARRAY[*]}; do
		wait $gzipid #wait for all gzip to finish
	done

	if $PAIRED_END; then
		cat $DOWNLOAD_DIR_TMP/*_1.fastq.gz > $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME"_1.fastq.gz"
		cat $DOWNLOAD_DIR_TMP/*_2.fastq.gz > $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME"_2.fastq.gz"
	else
		cat $DOWNLOAD_DIR_TMP/*.fastq.gz   > $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME".fastq.gz"
#	rm -r $DOWNLOAD_DIR_TMP
	fi
	printProgress "[masterDownload] Fastq files for $USER_INPUT_NAME are moved to $FASTQ_DIRECTORY"
}




# JULIEN: modify this to redirect the stdout to a log file
### Trim fastq files for adaptors using trimmomatic
### Could be called with $1=directory that holds fastq files
function trimReads () {	
	if [[ $TRIM_READ == false ]]; then #exit script if not needed, but make symbolic links to FASTQ files in the BAM folder
		printProgress "[trimReads] Read trimming not required. Skipping trim step."
		if [[ $PAIRED_END ]]; then
			checkFileExists $FASTQ_PATH"_1.fastq.gz"
			checkFileExists $FASTQ_PATH"_2.fastq.gz" 
			ln -s $FASTQ_PATH"_1.fastq.gz" $TEMP_DIR/$USER_INPUT_NAME"_1.fastq.gz"
			ln -s $FASTQ_PATH"_2.fastq.gz" $TEMP_DIR/$USER_INPUT_NAME"_2.fastq.gz"		
		else # single-end
			checkFileExists $FASTQ_PATH".fastq.gz" 
			ln -s $FASTQ_PATH".fastq.gz" $TEMP_DIR/$USER_INPUT_NAME".fastq.gz"
		fi
		return 0   # pass onto next data set
		
		
	elif [[ $TRIM_READ == true ]]; then
		cd $TEMP_DIR
		printProgress "[trimReads] Started at [$(date)]"

		if [[ $PAIRED_END == true ]]; then
			checkFileExists $FASTQ_PATH"_1.fastq.gz"
			checkFileExists $FASTQ_PATH"_2.fastq.gz"
			printProgress "[trimReads] Trimming "$USER_INPUT_NAME"_1.fastq.gz and "$USER_INPUT_NAME"_2.fastq.gz..."
			( $TRIMMOMATIC PE -threads $RUN_THREAD \
			-summary "$USER_INPUT_NAME"_trimSummary.txt \
			$FASTQ_PATH"_1.fastq.gz" $FASTQ_PATH"_2.fastq.gz" \
			$USER_INPUT_NAME"_1.fastq.gz" $USER_INPUT_NAME"_unpaired_trim_1.fastq.gz" $USER_INPUT_NAME"_2.fastq.gz" $USER_INPUT_NAME"_unpaired_trim_2.fastq.gz" \
			ILLUMINACLIP:$ILLUMINA_ADAPATORS_ALL":2:30:10" \
			SLIDINGWINDOW:4:20 \
			MINLEN:36 ) 2> "$USER_INPUT_NAME"_trimLog.txt
			cat "$USER_INPUT_NAME"_trimLog.txt "$USER_INPUT_NAME"_trimSummary.txt >> $LOG_FILE
	
			# Keep unpaired filtered reads for WGBS data
			if [[ $USER_INPUT_NAME == *"Bisulfite-Seq"* ]] || [[ $USER_INPUT_NAME == *"RRBS"* ]] || [[ $USER_INPUT_NAME == *"BSSeq"* ]] || [[ $USER_INPUT_NAME == *"PBAT"* ]] || [[ $USER_INPUT_NAME == *"DNAme"* ]]; then 
				return 0
			else
				rm $USER_INPUT_NAME"_unpaired_trim_1.fastq.gz" $USER_INPUT_NAME"_unpaired_trim_2.fastq.gz"
			fi

		else #Single-End
			checkFileExists $FASTQ_PATH".fastq.gz" 
			printProgress "[trimReads] Trimming $USER_INPUT_NAME.fastq.gz"
			( $TRIMMOMATIC SE -threads $RUN_THREAD \
			-summary "$USER_INPUT_NAME"_trimSummary.txt \
			$FASTQ_PATH".fastq.gz" $USER_INPUT_NAME".fastq.gz" \
			ILLUMINACLIP:$ILLUMINA_ADAPATORS_ALL":2:30:10" \
			SLIDINGWINDOW:4:20 \
			MINLEN:36 ) 2> "$USER_INPUT_NAME"_trimLog.txt
			cat "$USER_INPUT_NAME"_trimLog.txt "$USER_INPUT_NAME"_trimSummary.txt >> $LOG_FILE
		fi

		printProgress "[trimReads] Finished trimming fastq files for $USER_INPUT_NAME..."
		cd $CURRENT_DIRECTORY
	else
		printProgress "[trimReads] For some reason I don't know whether reads are PE or SE. Fix me!"
		
	fi

}


### Align fastq files to using data-specific aligner
### $1 can be a specific dataset prefix, will only align according to that tag
function masterAlign () {
	# if function iwithout an argument, use this argument instead
	SEARCH_KEY=${1:-$SEARCH_KEY}
	STAR_ARGUMENTS="--genomeDir $STAR_GENOME_DIR --runThreadN $RUN_THREAD --sjdbOverhang 70 --outFilterType BySJout --twopassMode Basic --twopass1readsN 1000000000 --outSAMunmapped Within --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --readFilesCommand zcat "
	BISMARK_ARGUMENTS="--chunkmbs $BISMARK_MEM -p $BOWTIE_THREAD --bowtie2 --bam $BISMARK_GENOME_DIR " 	# assumes that samtools and bowtie2 are in the path
	
	# -b
	if $BAM_INPUT; then #if input is BAM for trackhub - set up global variables and exit function
		# might have to set TEMP_DIR=$INPUT_BAM_DIRECTORY as for FASTQ_INPUT
		setupVariables $SEARCH_KEY
		return 0
	fi
	
	# -f
	if $FASTQ_INPUT; then
		TEMP_DIR=$FASTQ_DIRECTORY
	fi
	
	printProgress "[masterAlign] Starting..."

	# loop over replicates
	for FILE in $TEMP_DIR/*$SEARCH_KEY*fastq.gz; do

		setupVariables $FILE
		trimReads # this is skipped if set to false but nevertheless creates symbolic links of Fastq files into the BAM folder
		# I don't like this as some applications require custom read trimming (e.g. in PBAT you remove the first 4-8 Ns of each read)
		# maybe just trim the reads beforehand? In that case, teh function trimReads should take in $1 (SE) or $1 and $2 (PE)

		if [[ $FILE == *"RNA"* ]]; then
			alignSTAR
		elif [[ $FILE == *"RRBS"* ]] || [[ $FILE == *"BSSeq"* ]] || [[ $FILE == *"PBAT"* ]]; then
			alignBismark
		elif [[ $FILE == *"HiC" ]]; then
			alignHiCUP
		elif $ALLELE_SPECIFIC || $USE_BOWTIE; then
			alignBowtie2
		else
			alignBWA
		fi
	done
	printProgress "[masterAlign] Alignment of $SEARCH_KEY fastq files to $GENOME_BUILD completed."
	cd $CURRENT_DIRECTORY

}


function postAlignmentCleanup () {
	for BAM_FILE_TO_CLEAN in $TEMP_DIR/*$SEARCH_KEY*raw.bam; do
		printProgress "[postAlignmentCleanup] Started on $BAM_FILE_TO_CLEAN"
		if $ALLELE_RUN || [[ $BAM_FILE_TO_CLEAN == *"HiC" ]]; then #allele specific run will need to do unpacking
			continue #move onto next set instead of refining the produced RAW BAM
		fi	
		if [[ $BAM_FILE_TO_CLEAN == *"RRBS"* ]] || [[ $BAM_FILE_TO_CLEAN == *"BSSeq"* ]] || [[ $BAM_FILE_TO_CLEAN == *"PBAT"* ]]; then
			cleanAndExtractBismark $BAM_FILE_TO_CLEAN
		else # RNA-seq or ChIP-seq
			refineBam
		fi
	done
}



### Alignment of RNASeq data using STAR
function alignSTAR () {
	FILE_STAR_OUTPUT=$USER_INPUT_NAME"Aligned.out.bam"

	if $PAIRED_END; then
		printProgress "[masterAlign STAR] Aligning "$USER_INPUT_NAME"_1.fastq.gz and "$USER_INPUT_NAME"_2.fastq.gz to genome..."
		$STAR $STAR_ARGUMENTS --readFilesIn $FILE_FASTQ1 $FILE_FASTQ2 --outFileNamePrefix $USER_INPUT_NAME

	else #Single-End
		printProgress "[masterAlign STAR] Aligning $USER_INPUT_NAME.fastq.gz to genome..."
		$STAR $STAR_ARGUMENTS --readFilesIn $FILE_FASTQ --outFileNamePrefix $USER_INPUT_NAME
	fi
	checkFileExists $FILE_STAR_OUTPUT
	printProgress "[masterAlign STAR] Alignment completed -> $FILE_RAW_BAM"
	mv $FILE_STAR_OUTPUT $FILE_RAW_BAM
	rm -r $SEARCH_KEY*"STAR"*
}




### Alignment of HiC data with HiCUP, assumed to be paired end
function alignHiCUP () {
	ENZYME=$(grep -e $USER_INPUT_NAME $INPUT_SAMPLE_LIST_FILE | cut -f5 | head -n 1)

	if [[ -z $ENZYME ]]; then
		printProgress "[masterAlign HiCUP] ERROR:\tDigestion enzyme was not entered for $SEARCH_KEY HiC dataset..."
		exit 1
	fi

	printProgress "[masterAlign HiCUP] Digesting $GENOME_BUILD with $ENZYME"
	hicup_digester --re1 $ENZYME --outdir $TEMP_DIR	--genome $GENOME_BUILD"_"$SEARCH_KEY $GENOME_FILE
	local DIGEST_FILE=$(ls $TEMPDIR/Digest*$SEARCH_KEY*)

	printProgress "[masterAlign HiCUP] Aligning $USER_INPUT_NAME to $GENOME_BUILD"
	$HICUP --bowtie2 $BOWTIE2 --digest $DIGEST_FILE --index $GENOME_DIR/$GENOME_BUILD --outdir $TEMP_DIR --temp $TEMP_DIR --threads $RUN_THREAD $USER_INPUT_NAME"_1.fastq.gz" $USER_INPUT_NAME"_2.fastq.gz"

	printProgress "[masterAlign HiCUP] Alignment completed -> $FILE_RAW_BAM"
	rm $DIGEST_FILE
}

### Alignment of ChIPseq data with BWA-mem
function alignBWA () {
	FILE_SAM=$USER_INPUT_NAME".sam"

	if $PAIRED_END; then
		printProgress "[masterAlign BWA] Aligning "$USER_INPUT_NAME"_1.fastq.gz and "$USER_INPUT_NAME"_2.fastq.gz to genome..."
		$BWA mem -t $RUN_THREAD $GENOME_FILE $FILE_FASTQ1 $FILE_FASTQ2 > $FILE_SAM

	else #Single-End
		printProgress "[masterAlign BWA] Aligning "$USER_INPUT_NAME" to genome..."
		$BWA mem -t $RUN_THREAD $GENOME_FILE $FILE_FASTQ > $FILE_SAM
	fi

	printProgress "[masterAlign BWA] Converting sam file to .bam file..."
	$SAMTOOLS view -bhS -@ $RUN_THREAD $FILE_SAM > $FILE_RAW_BAM
	printProgress "[masterAlign BWA] Alignment completed -> $FILE_RAW_BAM"
	rm $FILE_SAM

}

### Alignment of ChIPSeq data using bowtie2 (for allelic specific)
function alignBowtie2 () {
	FILE_SAM=$USER_INPUT_NAME".sam"

	if $PAIRED_END; then
		printProgress "[masterAlign Bowtie2] Aligning "$USER_INPUT_NAME"_1.fastq.gz and "$USER_INPUT_NAME"_2.fastq.gz to genome..."
		$BOWTIE2 -x $BOWTIE2_INDEXES --local -p $RUN_THREAD -1 $FILE_FASTQ1 -2 $FILE_FASTQ2 -S $FILE_SAM

	else #Single-End
		printProgress "[masterAlign Bowtie2] Aligning $USER_INPUT_NAME.fastq to genome..."
		$BOWTIE2 -x $BOWTIE2_INDEXES --local -p $RUN_THREAD -U $FILE_FASTQ -S $FILE_SAM
	fi

	printProgress "[masterAlign Bowtie2] Converting sam file to .bam file..."
	$SAMTOOLS view -bhS -@ $RUN_THREAD $FILE_SAM > $FILE_RAW_BAM
	printProgress "[masterAlign Bowtie2] Alignment completed -> $FILE_RAW_BAM"
	rm $FILE_SAM

}

### Sort and mark duplicate alignments in BAM file
function refineBam () {
	local FILE_CLEANED_BAM=$USER_INPUT_NAME"_cleaned.bam"
	local FILE_SORTED_BAM=$USER_INPUT_NAME"_sort.bam"

	#soft clipping alignment that hangs off end of reference & set MAPQ to 0 for unmapped reads
	printProgress "[refineBAM] Refining $FILE_RAW_BAM"
	$PICARD CleanSam I=$FILE_RAW_BAM O=$FILE_CLEANED_BAM

	printProgress "[refineBAM] Sorting by coordinates..."
	$SAMTOOLS sort -@ $RUN_THREAD -m $THREAD_MEM -o $FILE_SORTED_BAM -T $USER_INPUT_NAME $FILE_CLEANED_BAM

	printProgress "[refineBAM] Marking duplicates..." #not removing the duplicates
	$PICARD MarkDuplicates I=$FILE_SORTED_BAM O=$FILE_BAM M=$USER_INPUT_NAME"_markDupeMetrics.txt"

	mv $FILE_BAM $OUTPUT_BAM_FOLDER
	printProgress "[refineBAM] Final $FILE_BAM is moved to $OUTPUT_BAM_FOLDER"

#	rm $FILE_RAW_BAM $FILE_CLEANED_BAM $FILE_SORTED_BAM #remove all the buffer bam files

}


### Searches for .bam files containing _Rep#.bam in nested directories and combines them into a single bam
function collapseReplicates () {
	cd $TEMP_DIR
	CURRENT=""
	printProgress "[collapseReplicates] Started at [$(date)]"
	
	# loop over replicates
	for FILE in $TEMP_DIR/*$SEARCH_KEY*.bam; do
		if [[ $FILE == *_[Rr]ep* ]]; then
			if [[ $CURRENT != ${FILE//_[Rr]ep*.bam/}*.bam ]]; then
				CURRENT=$FILE
				local REP_LIST=$(ls ${FILE//_[Rr]ep*.bam/}*.bam)
				local REPLICATE_COUNT=$(ls -l $REP_LIST | wc -l | awk '{print $1}' ) # ISNT THIS SIMPLER?
				local MERGED_BAM=${FILE//_[Rr]ep*.bam/_mergedReps1-"$REPLICATE_COUNT".bam}		
				
				printProgress "[collapseReplicates] Bam files to merge: $REP_LIST \n[collapseReplicates] Replicate count for $SEARCH_KEY: $REPLICATE_COUNT \n[collapseReplicates] Merged bam file name: $MERGED_BAM"
				$SAMTOOLS merge --threads $RUN_THREAD $MERGED_BAM $REP_LIST
				checkBamFileExists $MERGED_BAM
				printProgress "[collapseReplicates] Indexing and obtaining stats for replicates..."
				for x in $REP_LIST; do
					$SAMTOOLS index $x
					trueStats $x
				done
			fi

		else #file does not contain "rep or Rep"
			printProgress "[collapseReplicates] No replicates for $FILE"
			MERGED_BAM=$FILE
			$SAMTOOLS index $MERGED_BAM
			trueStats $MERGED_BAM
			checkBamFileExists $MERGED_BAM
		fi
	done

	printProgress "[collapseReplicates] All replicates for $SEARCH_KEY have been combined at [$(date)]"	
	cd $CURRENT_DIRECTORY
}


# to be run on duplicate-marked BAM files
# NOT compatible with DNAme data
function trueStats () {
	local INPUT_BAM=$1
	local TRUESTATS_OUTPUT=${INPUT_BAM//.bam/_trueStats.txt}
	checkFileExists $INPUT_BAM
	
	printProgress "[trueStats] Collecting alignment statistics for $INPUT_BAM"
	echo "" >> $TRUESTATS_OUTPUT
	echo "$INPUT_BAM trueStats" >> $TRUESTATS_OUTPUT
	echo "filtering $INPUT_BAM for mapped reads" >> $TRUESTATS_OUTPUT
	$SAMTOOLS view -b -F 4 $INPUT_BAM > $INPUT_BAM.mapped.bam
	$SAMTOOLS flagstat $INPUT_BAM.mapped.bam >> $TRUESTATS_OUTPUT
	echo "" >> $TRUESTATS_OUTPUT
	echo "filtering $INPUT_BAM for MAPQ >=$MIN_MAPQ" >> $TRUESTATS_OUTPUT
	$SAMTOOLS view -bq "$MIN_MAPQ" $INPUT_BAM.mapped.bam > $INPUT_BAM.mapped.MAPQ.bam
	$SAMTOOLS flagstat $INPUT_BAM.mapped.MAPQ.bam >> $TRUESTATS_OUTPUT
	rm $INPUT_BAM.mapped.bam $INPUT_BAM.mapped.MAPQ.bam
	cat $TRUESTATS_OUTPUT >> $LOG_FILE
}


### Checking dependencies of the functions
function checkDependencies () {
	printProgress "[checkDependencies] Checking Dependencies [$(date)]"
	EXIT=0
	for COMMAND in "${DEPENDENCIES[@]}"; do
		printProgress "[setup checkDependencies] $COMMAND..."
		command -v $COMMAND > /dev/null 2>&1 || {
			echo -e >&2 "\t\t$COMMAND not found!"
			EXIT=1
		}
	done

	if [[ $EXIT = 1 || $DEPEND = 1 ]] ; then
		exit 1
	fi
}


function removeFASTQ () {
	if [[ $KEEP_FASTQ == false ]]; then
		echo "Not keeping fastq..."
		local FASTQ_DIR=$CURRENT_DIRECTORY/$FASTQ_DIRECTORY
		rm $FASTQ_DIR/*$SEARCH_KEY*fastq.gz

		if [[ $(ls -1 $FASTQ_DIR | wc -l) == 0 ]]; then #if the Fastq directory isempty, then it can be removed
			rm -r $FASTQ_DIR
		fi

	fi
}


###Extract fastq file(s) from all bam within provided directory
function obtainFastqFromBAM () {
	cd $TEMP_DIR
	for BAM_INPUT in $CURRENT_DIRECTORY/${1:-$BAM_REALIGNMENT_DIRECTORY}/*.bam; do
		USER_INPUT_NAME=${BAM_INPUT##*/}
		SORTED="temp_sorted.bam"
		mkdir -p "temp"
		OUTPUT="temp"/${NAME//.bam/#.fastq} #the # will be replaced by _1/_2 for PE reads in bam2fastq

		echo "Sorting $USER_INPUT_NAME by read name..."
		$SAMTOOLS sort -@ $RUN_THREAD -m $THREAD_MEM -o $SORTED -n $BAM_INPUT
		echo "Extracting reads..."
		$BAM2FASTQ -o $OUTPUT $SORTED
		rm $SORTED

		echo "Compressing fastq file(s)..."
		for FILE in ${OUTPUT//#.fastq/*.fastq}; do
			( gzip $FILE
				mv $FILE.gz ${FILE//temp/$CURRENT_DIRECTORY\/$FASTQ_DIRECTORY}.gz
			) &
			COMPRESS_PID_ARRAY[$COMPRESS_COUNTER]=$!
			((COMPRESS_COUNTER++))
		done
	done
	cd $CURRENT_DIRECTORY
}


# move all of these to a separate .sh file?
########################################
#	DNA METHYLATION DATA FUNCTIONS     #
########################################

### Alignment of BSSeq data using Bismark
function alignBismark () {
	cd $TEMP_DIR
	BISMARK_OUTPUT=$USER_INPUT_NAME"_bismark_bt2.bam"
	if [[ $USER_INPUT_NAME == *"PBAT"* ]] ; then
		BISMARK_ARGUMENTS_LOCAL=$BISMARK_ARGUMENTS"--pbat "
	else
		BISMARK_ARGUMENTS_LOCAL=$BISMARK_ARGUMENTS"--non_directional " # uhhhhh are we sure this is correct? shouldn't non-pbat samples simply align with default?
	fi

	if [[ $PAIRED_END == true ]]; then
		# align PE reads
		printProgress "[masterAlign Bismark] Aligning "$USER_INPUT_NAME"_1.fastq.gz and "$USER_INPUT_NAME"_2.fastq.gz to genome..."
		# -un argument keeps mates that do not align, important for many WGBS libraries
		$BISMARK $BISMARK_ARGUMENTS_LOCAL -un -1 $TEMP_DIR/$USER_INPUT_NAME"_1.fastq.gz" -2 $TEMP_DIR/$USER_INPUT_NAME"_2.fastq.gz"
		BISMARK_OUTPUT=${BISMARK_OUTPUT//_bismark_bt2.bam/_1_bismark_bt2_pe.bam}
		mv $BISMARK_OUTPUT $TEMP_DIR/$USER_INPUT_NAME"_PE_raw.bam"
		checkBamFileExists $TEMP_DIR/$USER_INPUT_NAME"_PE_raw.bam" 
		
		# align unaligned reads too!
		if [[ $TRIM_READ == false ]]; then
			printProgress "[masterAlign Bismark] Aligning "$USER_INPUT_NAME" unaligned read 1 to genome..."
			$BISMARK $BISMARK_ARGUMENTS_LOCAL $TEMP_DIR/$USER_INPUT_NAME"_1.fastq.gz_unmapped_reads_1.fq.gz"
			printProgress "[masterAlign Bismark] Aligning "$USER_INPUT_NAME" unaligned read 2 to genome..."
			# do not run with --pbat option if it exists. Run default in order to align to all 4 reference genome strands
			$BISMARK $BISMARK_ARGUMENTS $TEMP_DIR/$USER_INPUT_NAME"_2.fastq.gz_unmapped_reads_2.fq.gz"
			rm $TEMP_DIR/$USER_INPUT_NAME"_1.fastq.gz_unmapped_reads_1.fq.gz" $TEMP_DIR/$USER_INPUT_NAME"_2.fastq.gz_unmapped_reads_2.fq.gz"

		else # read trimming occured, combine unpaired reads from trimming and from aligning
			printProgress "[masterAlign Bismark] Aligning "$USER_INPUT_NAME" unaligned read 1 and trimmed read 1 to genome..."
			cat $TEMP_DIR/$USER_INPUT_NAME"_unpaired_trim_1.fastq.gz" $TEMP_DIR/$USER_INPUT_NAME"_1.fastq.gz_unmapped_reads_1.fq.gz" > $TEMP_DIR/$USER_INPUT_NAME"_UNPAIRED_1.fastq.gz"
			$BISMARK $BISMARK_ARGUMENTS_LOCAL $TEMP_DIR/$USER_INPUT_NAME"_UNPAIRED_1.fastq.gz"
			rm $TEMP_DIR/$USER_INPUT_NAME"_unpaired_trim_1.fastq.gz" $TEMP_DIR/$USER_INPUT_NAME"_1.fastq.gz_unmapped_reads_1.fq.gz" $TEMP_DIR/$USER_INPUT_NAME"_UNPAIRED_1.fastq.gz"
			# outputs: $TEMP_DIR/$USER_INPUT_NAME"_UNPAIRED_1_bismark_bt2.bam"
			
			printProgress "[masterAlign Bismark] Aligning "$USER_INPUT_NAME" unaligned read 2 and trimmed read 2 to genome..."
			cat $TEMP_DIR/$USER_INPUT_NAME"_unpaired_trim_2.fastq.gz" $TEMP_DIR/$USER_INPUT_NAME"_2.fastq.gz_unmapped_reads_2.fq.gz" > $TEMP_DIR/$USER_INPUT_NAME"_UNPAIRED_2.fastq.gz"
			# do not run with --pbat option if it exists. Run default in order to align to all 4 reference genome strands
			$BISMARK $BISMARK_ARGUMENTS $TEMP_DIR/$USER_INPUT_NAME"_UNPAIRED_2.fastq.gz"
			rm $TEMP_DIR/$USER_INPUT_NAME"_unpaired_trim_2.fastq.gz" $TEMP_DIR/$USER_INPUT_NAME"_2.fastq.gz_unmapped_reads_2.fq.gz" $TEMP_DIR/$USER_INPUT_NAME"_UNPAIRED_2.fastq.gz"
			# outputs: $TEMP_DIR/$USER_INPUT_NAME"_UNPAIRED_2_bismark_bt2.bam"

		# combine SE bams
		printProgress "[masterAlign Bismark] Combining unpaired SE read Bams..."
		$SAMTOOLS merge --threads $RUN_THREAD $TEMP_DIR/$USER_INPUT_NAME"_SE_raw.bam" $TEMP_DIR/$USER_INPUT_NAME"_UNPAIRED_1_bismark_bt2.bam" $TEMP_DIR/$USER_INPUT_NAME"_UNPAIRED_2_bismark_bt2.bam"
		rm $TEMP_DIR/$USER_INPUT_NAME"_UNPAIRED_1_bismark_bt2.bam" $TEMP_DIR/$USER_INPUT_NAME"_UNPAIRED_2_bismark_bt2.bam"
		fi
		
	else #Single-End
		printProgress "[masterAlign Bismark] Aligning $USER_INPUT_NAME.fastq.gz to genome..."
		$BISMARK $BISMARK_ARGUMENTS_LOCAL $TEMP_DIR/$USER_INPUT_NAME".fastq.gz"
		mv $BISMARK_OUTPUT $TEMP_DIR/$USER_INPUT_NAME"_SE_raw.bam"
		checkBamFileExists $TEMP_DIR/$USER_INPUT_NAME"_SE_raw.bam" 
	fi
	
	BISMARK_ARGUMENTS_LOCAL="" # reset variable to avoid compounding it (e.g. --pbat --pbat --pbat in the case of 3 replicate samples)
	
	echo "Bismark alignment statistics"  >> $LOG_FILE
	cat $TEMP_DIR/$USER_INPUT_NAME*bt2*report.txt >> $LOG_FILE
	printProgress "[masterAlign Bismark] Alignment of $USER_INPUT_NAME completed!"
}




# EXTRACT STEP TAKES A LONG TIME. DO NOT CODE MORE STUFF INTO THIS FUNCTION
function cleanAndExtractBismark () {
	# this function is run on all files ending in *raw.bam
	# NOTE TO JRA (HIMSELF): BISMARK DEDUPLICATE IS EQUIVALENT TO PICARD MARKDUPLICATES PLUS SAMTOOLS FILTERING (WITH -F 1024)
	local CLEANBIS_INPUT=$1
	cd $TEMP_DIR
	
	printProgress "[cleanAndExtractBismark] Deduplicating reads..."
	if [[ $CLEANBIS_INPUT == *_PE_* ]]; then
		if [[ $CLEANBIS_INPUT == *RRBS* ]]; then
			DEDUPLICATED_BAM=$CLEANBIS_INPUT # it is not recommended to deduplicate reads generated by RRBS
		else
			DEDUPLICATED_BAM=${CLEANBIS_INPUT//.bam/.deduplicated.bam}
			$BISMARK_DEDUPLICATE -p $CLEANBIS_INPUT --bam # outputs file ${FILE_RAW_BAM//.bam/.deduplicated.bam} and ${FILE_RAW_BAM//.bam/.deduplicated_report.txt}
			cat ${DEDUPLICATED_BAM//deduplicated.bam/deduplication_report.txt} >> $LOG_FILE
		fi
		$BISMARK_METH_EXTRACT -p --multicore $RUN_THREAD --comprehensive --merge_non_CpG --bedGraph --counts --buffer_size $THREAD_MEM --cytosine_report --genome_folder $GENOME_DIR $DEDUPLICATED_BAM
		#outputs .CpG_report.txt
		
	elif [[ $CLEANBIS_INPUT == *_SE_* ]]; then
		if [[ $CLEANBIS_INPUT == *RRBS* ]]; then
			DEDUPLICATED_BAM=$CLEANBIS_INPUT
		else
			DEDUPLICATED_BAM=${CLEANBIS_INPUT//.bam/.deduplicated.bam}
			$BISMARK_DEDUPLICATE -s $CLEANBIS_INPUT --bam
			cat ${DEDUPLICATED_BAM//deduplicated.bam/deduplication_report.txt} >> $LOG_FILE
		fi
		$BISMARK_METH_EXTRACT -s --multicore $RUN_THREAD --comprehensive --merge_non_CpG --bedGraph --counts --buffer_size $THREAD_MEM --cytosine_report --genome_folder $GENOME_DIR $DEDUPLICATED_BAM		

	else 
		printProgress "[cleanAndExtractBismark] Unknown file type (name should contain SE or PE). Exiting."
		exit 0

	fi
	
	# COLLECT SOME BASIC STATS
	#collectBismarkDuplicationCoverageAndConversionStats $DEDUPLICATED_BAM
	
}


### merge a cytosine report into a CpG site report
function mergeTwoStrandMethylation () {
	local CPG_REPORT_FILE=$1
	checkFileExists $CPG_REPORT_FILE
	printProgress "[mergeTwoStrandMethylation] Started mergeTwoStrandMethylation on $CPG_REPORT_FILE"

	sort -k1,1 -k2,2n $CPG_REPORT_FILE > ${CPG_REPORT_FILE//.txt/_sort.txt}
	awk '
	BEGIN {
			FS = "\t"
			OFS = "\t"

			CHR = ""
			FIRST_POS = 0
			METHYL = 0
			UNMETHYL = 0
			FIRST_TRI = ""
			START = 0
	}
	{
		if ($2 == FIRST_POS || $2 == FIRST_POS + 1) {
			METHYL += $4
			UNMETHYL += $5
		} else {
			if (START != 0 ) {
				if (METHYL + UNMETHYL > 0) {
					printf CHR "\t" FIRST_POS "\t" 
					printf "%6f\t", METHYL / (METHYL + UNMETHYL) * 100.0
					print METHYL, UNMETHYL, $6, FIRST_TRI
				} else {
					print CHR, FIRST_POS, "NA", METHYL, UNMETHYL, $6, FIRST_TRI
				}
			}
			START = 1
			CHR = $1
			FIRST_POS = $2
			METHYL = $4
			UNMETHYL = $5
			FIRST_TRI = $7
		}
	}
	END {
		if (METHYL + UNMETHYL > 0) {
			printf CHR "\t" FIRST_POS "\t" 
			printf "%6f\t", METHYL / (METHYL + UNMETHYL) * 100.0
			print METHYL, UNMETHYL, $6, FIRST_TRI
		} else {
			print CHR, FIRST_POS, "NA", METHYL, UNMETHYL, $6, FIRST_TRI
		}
	}' ${CPG_REPORT_FILE//.txt/_sort.txt} > ${CPG_REPORT_FILE//.txt/_mergeTwoStrands.txt}
	rm ${CPG_REPORT_FILE//.txt/_sort.txt}
	checkFileExists ${CPG_REPORT_FILE//.txt/_mergeTwoStrands.txt}
}


# If _PE_ and _SE_ exist, combine and rename to _SEPE_, and remove uncombined files
function collapseSEPEbismark () {
	
	cd $TEMP_DIR
	printProgress "[collapseSEPEbismark] Started at [$(date)]"
	for FILE in $TEMP_DIR/*$SEARCH_KEY*.CpG_report.txt; do
		# Combine SE and PE if PE exists
		if [[ $FILE == *_PE_* ]]; then
			SE_FILE=${FILE//_PE_/_SE_}
			# Check if the unpaired SE file exists (PBAT)
			if [[ -f $SE_FILE ]]; then
				printProgress "[collapseSEPEbismark] Combining $FILE + $SE_FILE -> ${FILE//_PE_/_SEPE_}"
				cat $FILE $SE_FILE > ${FILE//_PE_/_SEPE_}
				rm  $FILE $SE_FILE
			else # do not do anything if file is _PE_ but --un was not set duringh alignment (not _SE_)
				continue
			fi
		fi # do nothing if file does not contain _PE_
	done
	
	# get TrueStats because I love TrueStats. Then remove all intermediate files created by this loop
	for FILE in $TEMP_DIR/*$SEARCH_KEY*raw.bam; do
		if [[ $FILE == *_PE_* ]]; then
			SE_FILE=${FILE//_PE_/_SE_}
			if [[ -f $SE_FILE ]]; then
				printProgress "[collapseSEPEbismark] Combining $FILE + $SE_FILE -> ${FILE//_PE_/_SEPE_}"
				$SAMTOOLS merge --threads $RUN_THREAD ${FILE//_PE_/_SEPE_} $FILE $SE_FILE
				mv $FILE $OUTPUT_BAM_FOLDER/
				mv $SE_FILE $OUTPUT_BAM_FOLDER/				
			else
				continue
			fi
		fi
	done

	for FILE in $TEMP_DIR/CpG*$SEARCH_KEY*; do
		if [[ $FILE == *_PE_* ]]; then
			SE_FILE=${FILE//_PE_/_SE_}
			if [[ -f $SE_FILE ]]; then
				printProgress "[collapseSEPEbismark] Combining $FILE + $SE_FILE -> ${FILE//_PE_/_SEPE_}"
				cat $FILE $SE_FILE  > ${FILE//_PE_/_SEPE_}		
				local NON_CPG_REPORT=${FILE//CpG_context_/Non_CpG_context_}
				printProgress "[collapseSEPEbismark] Combining $NON_CPG_REPORT + ${NON_CPG_REPORT//_PE_/_SE_} -> ${FILE//_PE_/_SEPE_}"
				cat $NON_CPG_REPORT ${NON_CPG_REPORT//_PE_/_SE_} > ${NON_CPG_REPORT//_PE_/_SEPE_} 			
			else
				continue
			fi
		fi
	done

	printProgress "[collapseSEPEbismark] Done at [$(date)]"
	
}


function collectBismarkStats () {
	cd $TEMP_DIR
	
	printProgress "[collectBismarkStats] Collecting bisulphite conversion rate and coverage for files $(ls *raw.bam)"
	
	for FILE in $TEMP_DIR/$SEARCH_KEY*raw.bam; do
		$PICARD SortSam I=$FILE O=${FILE//.bam/_tmp1.bam} SO=coordinate
		rm $FILE
		$PICARD MarkDuplicates I=${FILE//.bam/_tmp1.bam} O=$FILE M=trash.txt
		trueStats $FILE

		$SAMTOOLS view -bh -F 1540 $FILE > ${FILE//.bam/.rmDupbam}
		$SAMTOOLS coverage ${FILE//.bam/.rmDupbam} | grep -v "random" | grep -v "Un" >> $LOG_FILE
		$SAMTOOLS coverage ${FILE//.bam/.rmDupbam} | grep -v "random" | grep -v "Un" | grep -e "chr[1-9]" | awk '{sum+=$7} END { print "Autosomal avg coverage = ",sum/NR}' >> $LOG_FILE
		rm $FILE ${FILE//.bam/_tmp1.bam} ${FILE//.bam/.rmDupbam} trash.txt
		
	done

	for FILE in $TEMP_DIR/CpG*$SEARCH_KEY*; do
		local NON_CPG_REPORT=${FILE//CpG_context_/Non_CpG_context_}
				
		printProgress "[collapseSEPEbismark] Counting the number of lamba methylated (Z,X,H) and unmethylated (z,x,h) CpGs in $FILE and $NON_CPG_REPORT"
		echo "Conversion rate(%) = converted Cs/(converted+unconverted Cs)" >> ${FILE//CpG/conversion_tmp}
		cat $FILE $NON_CPG_REPORT | grep -e "J02459.1" | cut -f5 | sort | uniq -c >> ${FILE//CpG/conversion_tmp}
		awk '{ if ($2=="Z" || $2=="H" || $2=="X") summet+=$1; else sumunmet+=$1} END {print "Conversion rate: "sumunmet/(sumunmet+summet)"\tUnconverted (methylated): "summet"\tConverted (unmethylated): "sumunmet}' ${FILE//CpG/conversion_tmp} >> $LOG_FILE
		cat ${FILE//CpG/conversion_tmp} >> $LOG_FILE
		rm  $FILE $NON_CPG_REPORT ${FILE//CpG/conversion_tmp}
		
	done
}

function collapseReplicatesBismark () {
	cd $TEMP_DIR
	printProgress "[collapseReplicatesBismark] Started at [$(date)]"
	CURRENT=""
	local SUFFIX=".CpG_report.txt"
	for FILE in $TEMP_DIR/*$SEARCH_KEY*$SUFFIX; do
		if [[ $FILE == *_[Rr]ep* ]]; then
			if [[ $CURRENT != ${FILE//_[Rr]ep*$SUFFIX/}*$SUFFIX ]]; then
				CURRENT=$FILE
				local REP_LIST=$( ls ${FILE//_[Rr]ep*$SUFFIX/}*$SUFFIX )
				local REPLICATE_COUNT=$( ls -l $REP_LIST | wc -l | awk '{print $1}' )
				local MERGED_CPG_REPORT_TMP=$( echo $REP_LIST | cut -d ' ' -f1 )
				local MERGED_CPG_REPORT_TMPE=${MERGED_CPG_REPORT_TMP//_[Rr]ep[123456789]_/_} # replace repX with reps1-N
				local MERGED_CPG_REPORT=${MERGED_CPG_REPORT_TMPE//$SUFFIX/_mergedReps1-"$REPLICATE_COUNT"$SUFFIX}
				printProgress "[collapseReplicatesBismark] CpG_report.txt files to merge: $REP_LIST \n[collapseReplicatesBismark] Replicate count for $SEARCH_KEY: $REPLICATE_COUNT \n[collapseReplicatesBismark] Merged bam file name: $MERGED_CPG_REPORT"
				printProgress "[collapseReplicatesBismark] Merging and sorting CpG_report replicates..."
				# sorting takes about 3 minutes per CpG report on mbp16
				cat $REP_LIST > $MERGED_CPG_REPORT
				mergeTwoStrandMethylation $MERGED_CPG_REPORT
				for x in $REP_LIST; do
					mergeTwoStrandMethylation $x
					rm $x
				done
			fi
		else #file does not contain "rep or Rep"
			printProgress "[collapseReplicatesBismark] No replicates for $FILE"
		fi
	done
	printProgress "[collapseReplicatesBismark] Done at [$(date)]"
	
	
	
}	


function convertMethylationToBigWig () {
	local INPUT_CPG_REPORT=$1
	local VAR_MIN_DEPTH=$2
	printProgress "[masterTrackHub] convertMethylationToBigWig $INPUT_CPG_REPORT x$VAR_MIN_DEPTH"

	awk -v MIN_DEPTH=$VAR_MIN_DEPTH '
		BEGIN {
			print "#chr\tstart\tend\tstrand\tvalue"
			FS = "\t"
			OFS = "\t"
			if (MIN_DEPTH < 1)
				MIN_DEPTH = 1
		}
	$5 + $6 >= MIN_DEPTH {
			print $1, $2, $2+1, $3
	}' $INPUT_CPG_REPORT >  ${INPUT_CPG_REPORT/.txt/_x$VAR_MIN_DEPTH.bedGraph}

	$BEDGRAPHTOBW ${INPUT_CPG_REPORT/.txt/_x$VAR_MIN_DEPTH.bedGraph} $CHROM_SIZES ${INPUT_CPG_REPORT/.txt/_x$VAR_MIN_DEPTH.bw}
	# count number of covered CpGs
	CPG_COUNT=$( grep -v "J02459.1" ${INPUT_CPG_REPORT/.txt/_x$VAR_MIN_DEPTH.bedGraph}  |  wc -l  |  awk '{print $1}' )
	echo "$CPG_COUNT CpGs covered by $VAR_MIN_DEPTH reads in $INPUT_CPG_REPORT" >> $LOG_FILE
	rm ${INPUT_CPG_REPORT/.txt/_x$VAR_MIN_DEPTH.bedGraph}
	
}



function cleanupBismark  () {
	printProgress "[cleanupBismark] Moving relevant output files to $CURRENT_DIRECTORY"
	mv $TEMP_DIR/$SEARCH_KEY*.bw $OUTPUT_BIGWIG_FOLDER/
	mv $TEMP_DIR/$SEARCH_KEY*CpG_report_mergeTwoStrands.txt $OUTPUT_CPG_REPORT_FOLDER/
	cp $LOG_FILE $OUTPUT_STATS_FOLDER/
	
}




#########################################
#		TRACK-HUB FUNCTIONS		     	#
#########################################

function masterTrackHub () {
	
	BAM_COVERAGE_ARGUMENTS="--binSize $BIN_SIZE -p $RUN_THREAD --normalizeUsing $NORMALIZE --smoothLength $SMOOTH_WIN --outFileFormat bigwig --minMappingQuality $MIN_MAPQ --ignoreDuplicates"
	PRINTED_DIR=""
	# special little condition for DNAme data
	if [[ $USER_INPUT_NAME == *"BSSeq"* ]] || [[ $USER_INPUT_NAME == *"RRBS"* ]] || [[ $USER_INPUT_NAME == *"PBAT"* ]]; then
		
		for CPG_REPORT_FILE in $TEMP_DIR/*$SEARCH_KEY*CpG_report_mergeTwoStrands.txt; do
			convertMethylationToBigWig $CPG_REPORT_FILE 1
			convertMethylationToBigWig $CPG_REPORT_FILE 2
			convertMethylationToBigWig $CPG_REPORT_FILE 3
			convertMethylationToBigWig $CPG_REPORT_FILE 4
			convertMethylationToBigWig $CPG_REPORT_FILE 5
			convertMethylationToBigWig $CPG_REPORT_FILE 10
			convertMethylationToBigWig $CPG_REPORT_FILE 15
			convertMethylationToBigWig $CPG_REPORT_FILE 20
		done


	# not DNAme
	# RECODE THIS
	else

		for BAM_FILE in $TEMP_DIR/*$SEARCH_KEY*.bam; do
#		for BAM_FILE in ./*/*$SEARCH_KEY*.bam; do #currently in $CURRENT_DIRECTORY
			FOLDER_FILE=${BAM_FILE//.\//} #getting rid of the "./"
			FILE=$(basename $BAM_FILE) #leaving just the basename of the file
			FILE_NAME=${FILE//.bam/_}$NORMALIZE ##basename without file extension

			if [[ $FILE == *"BSSeq"* ]] || [[ $FILE == *"RRBS"* ]] || [[ $FILE == *"PBAT"* ]]; then
				break
			fi


			if [[ $BIN_SIZE != 1 ]] ; then
				FILE_NAME=$FILE_NAME"_b"$BIN_SIZE
			fi

			if [[ $SMOOTH_WIN != 0 ]] ; then
				FILE_NAME=$FILE_NAME"_s"$SMOOTH_WIN
			fi

			FOLDER_NAME=${FOLDER_FILE%%\/*} #removing longest text of the matching pattern (after "/" in this case)

			if [[ $FOLDER_NAME != $PRINTED_DIR ]] ; then # Create supertrack for housing associated tracks
				printf "track $FOLDER_NAME \nsuperTrack on show\nshortLabel $FOLDER_NAME \nlongLabel $FOLDER_NAME \n\n" | tee -a $TRACKDB
				PRINTED_DIR=$FOLDER_NAME
			fi

			FLAG=$($SAMTOOLS view $FOLDER_FILE | head -n 1 | cut -f2)
			if [[ $((FLAG&1)) == 1 ]] ; then #Paired-end
				PAIRED=true
			else
				PAIRED=false
			fi
			
			echo "Data are Paired-end:" $PAIRED

			if [[ $FILE == *"RNA"* ]]; then
				generateRNATrack	
			#ChIPseq - not stranded
			else
				generateBigwigsUnstranded $FOLDER_FILE $FILE_NAME
				printTrackHubUnstranded $FOLDER_NAME $FILE_NAME
			fi
			
		rm Actb.bed
		done
	fi

}


function generateRNATrack () {
	if [[ $PAIRED == true ]] ; then
		echo "Extracting F reads over Actb..."
		$SAMTOOLS view -L Actb.bed -f 64 $FOLDER_FILE > Actb.sam
	else
		echo "Extracting reads over Actb..."
		$SAMTOOLS view -L Actb.bed $FOLDER_FILE > Actb.sam
	fi

	STRANDED=$(awk 'BEGIN{PLUS=0; MINUS=0} {
		if( and($2,16) == 16) {
					MINUS++;
				} else {
					PLUS++;
				}
			} END{
				if( MINUS/(MINUS+PLUS) > 0.9) {
					print "Same-Strand";
				} else {
					if( MINUS/(MINUS+PLUS) < 0.1) {
						print "Opposite-Strand";
					} else {
						print "Unstranded";
					}
				}
			}' Actb.sam)
	rm Actb.sam
	echo "Data are" $STRANDED

	if [[ $STRANDED == "Unstranded" ]] ; then
		generateBigwigsUnstranded $FOLDER_FILE $FILE_NAME
		printTrackHubUnstranded $FOLDER_NAME $FILE_NAME
	else
		FILE_BIGWIG_POS=$TRACK_FOLDER"$FILE_NAME""_pos.bw"
		FILE_BIGWIG_NEG=$TRACK_FOLDER"$FILE_NAME""_neg.bw"
		FILE_BIGWIG_TEMP=$TEMP_DIR"/temp.bw"
		FILE_BEDGRAPH_TEMP=$TEMP_DIR"/temp.bedgraph"
		FILE_BEDGRAPH_TEMP2=$TEMP_DIR"/temp2.bedgraph"
		$BAMCOVERAGE $BAM_COVERAGE_ARGUMENTS -b $FOLDER_FILE --filterRNAstrand forward --outFileName $FILE_BIGWIG_POS
		$BAMCOVERAGE $BAM_COVERAGE_ARGUMENTS -b $FOLDER_FILE --filterRNAstrand reverse --outFileName $FILE_BIGWIG_NEG
		if [[ $STRANDED == "Opposite-Strand" ]] ; then
			mv $FILE_BIGWIG_NEG $FILE_BIGWIG_TEMP
			mv $FILE_BIGWIG_POS $FILE_BIGWIG_NEG
			mv $FILE_BIGWIG_TEMP $FILE_BIGWIG_POS
		fi
		bigWigToBedGraph $FILE_BIGWIG_NEG $FILE_BEDGRAPH_TEMP
		awk -F "\t" 'BEGIN {OFS="\t" } {
			$4 = 0 - $4;
			print $1, $2, $3, $4;
		}' $FILE_BEDGRAPH_TEMP > $FILE_BEDGRAPH_TEMP2
		mv $FILE_BEDGRAPH_TEMP2 $FILE_BEDGRAPH_TEMP
		$BEDGRAPHTOBW $FILE_BEDGRAPH_TEMP $CHROM_SIZES $FILE_BIGWIG_NEG
		rm $FILE_BEDGRAPH_TEMP
		printTrackHubStranded $FOLDER_NAME $FILE_NAME
	fi
}


#$1=Name of filtered .bam file $2=Final name
function generateBigwigsUnstranded () {
	echo "Generating bigwig files..."
	$SAMTOOLS index $1
	$BAMCOVERAGE $BAM_COVERAGE_ARGUMENTS -b $1 --outFileName $TRACK_FOLDER"$2.bw"

}

#$1=Supertrack name $2=Track name
function printTrackHubUnstranded () {
	printf "\ttrack %s\n\tparent %s\n\tshortLabel %s\n\tlongLabel %s\n\ttype bigWig\n\tbigDataUrl %s\n\tcolor 200,50,0\n\tvisibility full\n\tmaxHeightPixels 100:60:25\n\tautoScale on\n\talwaysZero on\n\n" $2 $1 $2 $2 $2".bw" | tee -a $TRACKDB

}

#Takes in supertrack name then track name
function printTrackHubStranded () {
	printf "\ttrack %s\n\tparent %s\n\tcontainer multiWig\n\tshortLabel %s\n\tlongLabel %s\n\ttype bigWig\n\tvisibility full\n\tmaxHeightPixels 100:60:25\n\tconfigurable on\n\tautoScale on\n\talwaysZero on\n\taggregate transparentOverlay\n\tshowSubtrackColorOnUi on\n\tpriority 1.0\n\n" $2 $1 $2 $2 | tee -a $TRACKDB
	printf "\t\ttrack %s\n\t\tparent %s\n\t\tshortLabel %s\n\t\tlongLabel %s\n\t\ttype bigWig\n\t\tbigDataUrl %s\n\t\tcolor 200,50,0\n\t\tautoScale on\n\n" $2"_pos" $2 $2"_pos" $2"_pos" $2"_pos.bw" | tee -a $TRACKDB
	printf "\t\ttrack %s\n\t\tparent %s\n\t\tshortLabel %s\n\t\tlongLabel %s\n\t\ttype bigWig\n\t\tbigDataUrl %s\n\t\tcolor 200,50,0\n\t\tautoScale on\n\n" $2"_neg" $2 $2"_neg" $2"_neg" $2"_neg.bw" | tee -a $TRACKDB

}
	
	
	
	
	
	
	
	
########################################
#	ALLELE-SPECIFIC FUNCTIONS      #
########################################

function checkPseudogenome () {
	if $ALLELE_SPECIFIC && $SEP_PARA ; then

		cd $TEMP_DIR

		CROSS_LIST=$(awk '(NF==4 && ($3!="" || $4!="")) {print $3"_"$4}' $INPUT_SAMPLE_LIST_FILE | sort | uniq) #create array of crosses
		EXIT_SCRIPT=false

		for CROSS in $CROSS_LIST; do
			HAPLO_1=${CROSS//_*}
			HAPLO_2=${CROSS##*_}

			if [[ -z $HAPLO_1 || -z $HAPLO_2 ]]; then
				echo -e "ERROR:\t One of your samples only has 1 haplotype entered. \nPlease to enter both haplotypes." \
				| tee -a
				EXIT_SCRIPT=true
				continue #continue to next iteration of cross
			fi

			MAKE_GENOME=true

			for DIP in $DIPLOID_GENOME_DIR/*; do #searching for folder with name of both haplotypes
				if [[ $DIP == *$HAPLO_1* && $DIP == *$HAPLO_2* ]]; then
					MAKE_GENOME=false

					if [[ ! -f $HAPLOID_GENOME_DIR/$HAPLO_1/$HAPLO_1".fa.refmap" ]] || \
						 [[ ! -f $HAPLOID_GENOME_DIR/$HAPLO_2/$HAPLO_2".fa.refmap" ]]; then #check RefMaps exist
						echo -e "ERROR:\t Either refmap file for haplotype $HAPLO_1 or $HAPLO_2 doesn't exists."
						EXIT_SCRIPT=true
					fi

					echo "Reference genome for $CROSS is found"
					break #exit the loop of iterating thru all the diploid genomes for this combo of cross
				fi
			done

			if $MAKE_GENOME; then
				echo "ERROR:\t Diploid genome was not found for $CROSS. \nPlease use provided CreateDipPseudoGenome.sh to obtain diploid pseudogenome."
				EXIT_SCRIPT=true
			fi

		done

		cd $CURRENT_DIRECTORY

		if $EXIT_SCRIPT; then
			echo -e "Errors were found with setup. \nPlease check error message before running again."
			exit #exit script to ensure input is formatted correctly before running the whole script
		fi
	fi
}

function setPseudogenome () {
	#NO!
	USER_INPUT_NAME=$SEARCH_KEY
	printProgress "[setPseudogenome] Started at [$(date)]"

	HAPLO_1=$(grep -e $USER_INPUT_NAME $INPUT_SAMPLE_LIST_FILE | cut -f3 | head -n 1)
	HAPLO_2=$(grep -e $USER_INPUT_NAME $INPUT_SAMPLE_LIST_FILE | cut -f4 | head -n 1)

	if [[ -z $HAPLO_1 || -z $HAPLO_2 ]]; then #if either haplotype is not entered, no allele specific for this data
		printProgress "[setPseudogenome] Haplotypes not provided, will not run allele specific pipeline on $USER_INPUT_NAME"
		continue #dont continue the function for this SRACODE entry (move on to next code)
	fi

	CROSS_COMBO=$HAPLO_1"_"$HAPLO_2
	printProgress "[setPseudogenome] Haplotypes for $USER_INPUT_NAME are $HAPLO_1 and $HAPLO_2."

	for DIP in $DIPLOID_GENOME_DIR/*; do #searching for folder with name of both haplotypes
		if [[ $DIP == *$HAPLO_1* && $DIP == *$HAPLO_2* ]]; then
			GENOME_BUILD=$(basename $DIP)
			GENOME_FILE=$DIP/*".fa"
			STAR_GENOME_DIR=$DIP/*"-STAR"
			BOWTIE2_INDEXES=$DIP/$GENOME_BUILD
			BISMARK_INDEXES=$DIP/Bisulfite_Genome #added by Julien
			HAPLO_1_REFMAP=$HAPLOID_GENOME_DIR/$HAPLO_1/$HAPLO_1".fa.refmap"
			HAPLO_2_REFMAP=$HAPLOID_GENOME_DIR/$HAPLO_2/$HAPLO_2".fa.refmap"
			printProgress "[setPseudogenome] Reference directory for $USER_INPUT_NAME is set to $DIP"
			break
		fi
	done
}

function unpackAllelic () { #working on bam that has already aligned to the pseudogenome (once per haplotype)
	cd $TEMP_DIR
	local HAPLO=$1

	printProgress "[unpackAllelic] Unpacking reads aligned to $HAPLO with MAPQ = $MIN_MAPQ started at [$(date)]"

	for TOT_RAW_BAM in $TEMP_DIR/$SEARCH_KEY*"_raw.bam"; do #should be in replicates (NAME1_rep1_raw.bam)
	#take SAM file, align to pseudogenome
	#split header into two, to separate the reads into two haplotypes, rename chr from hap1_chr to chr
		NAME_HAPLO_SUFFIX=$SEARCH_KEY"_"$HAPLO"_q$MIN_MAPQ"${TOT_RAW_BAM//$SEARCH_KEY/} #inserting haplotype name (-> NAME1_HAPLO_qMINMAPQ_rep1_raw.bam)
		#NONONO
		USER_INPUT_NAME=${NAME_HAPLO//_raw.bam/} #will include _rep if applicable
		FILE_RAW_BAM=$USER_INPUT_NAME"_raw.bam"

		printProgress "[unpackAllelic $HAPLO] Obtaining haplotype-specific header..."
		#take lines that include the haplotype name - haplotype_specific chrom & commands
		$SAMTOOLS view -H $TOT_RAW_BAM \
		| awk '( $0 ~ "'$HAPLO'" ) {print $0}' \
		| sed 's/'$HAPLO'_chr/chr/g' > $FILE_RAW_BAM

		# get UNIQUELY ALIGNED READS, separate into two files, only keep reads where their mate also maps to the same chromo of the same haplotype
		# uniquely aligned defined by MIN_MAPQ
		printProgress "[unpackAllelic $HAPLO] Obtaining haplotype-specific reads with MAPQ >= $MIN_MAPQ"
		$SAMTOOLS view $TOT_RAW_BAM -q $MIN_MAPQ \
		| awk '(( $3 ~ "'$HAPLO'" ) && ( $7 ~ "'$HAPLO'" || $7 == "*" || $7 == "=" )) {print $0}' \
		| sed 's/'$HAPLO'_chr/chr/g' >> $FILE_RAW_BAM

		printProgress "[unpackAllelic $HAPLO] Finished unpacking $HAPLO specific reads from $TOT_RAW_BAM -> $FILE_RAW_BAM"

		refineBAM # does this not take an argument here?

#		rm $TOT_RAW_BAM # need to keep this in order to extract reads on HAPLO-2
	done

	printProgress "[unpackAllelic] Unpacking bams for $HAPLO completed at [$(date)]"
	cd $CURRENT_DIRECTORY
}

function projectAllelic () {
	if [[ $ALLELE_SPECIFIC == false ]]; then
		return 0 #dont run this code if the script is not allele specific
	fi

	cd $TEMP_DIR

	printProgress "[projectAllelic] Started at [$(date)]"

	for FILE_BAM in $TEMP_DIR/*".bam"; do
		if [[ $FILE_BAM != *$HAPLO_1* || $FILE_BAM != *$HAPLO_2* ]]; then #only projecting bams that were aligned to pseudogenome
			continue
		fi

		printProgress "[projectAllelic] Removing duplicates from $FILE_BAM with F=$FLAG"
		local NAME_MAP_FLAG=${FILE_BAM//.bam/}"_F"$FLAG
		FLAG_BAM=$USER_INPUT_NAME_MAP_FLAG".bam"
		$SAMTOOLS view -bh -F $FLAG $FILE_BAM > $FLAG_BAM

		#nonscaled projection
		printProgress "[projectAllelic] Converting $FLAG_BAM to bedGraph"
		$BEDTOOLS genomecov -ibam $FLAG_BAM -bg -split -scale $SCALING_FACTOR > $USER_INPUT_NAME_MAP_FLAG"_preProject.bedgraph"
		prepWigAndProject $USER_INPUT_NAME_MAP_FLAG $USER_INPUT_NAME_MAP_FLAG"_preProject.bedgraph"

		#RPM scaled projection
		local READ_COUNT=$(samtools view -c $FLAG_NAME".bam")
		local SCALING_FACTOR=$(echo "scale=25; 1000000/$READ_COUNT" | bc) #calculating with 25 decimal places at least
		printProgress "[projectAllelic RPM] Detected $READ_COUNT filtered reads"

		if [[ $STRANDED_ALLELIC ]]; then #stranded RPM
			printProgress "[projectAllelic stranded RPM] Splitting $FLAG_BAM reads by 1st in pair..."
			local FIRST_PAIR_BAM=$USER_INPUT_NAME_MAP_FLAG"_firstInPair.bam"
			$SAMTOOLS view -bh -f 0x0040 $FLAG_BAM > $FIRST_PAIR_BAM
			$BEDTOOLS genomecov -ibam $FIRST_PAIR_BAM -bg -split -scale $SCALING_FACTOR -strand + > $USER_INPUT_NAME_MAP_FLAG"_first_pos_RPM.bedGraph"
			$BEDTOOLS genomecov -ibam $FIRST_PAIR_BAM -bg -split -scale $SCALING_FACTOR -strand - > $USER_INPUT_NAME_MAP_FLAG"_first_neg_RPM.bedGraph"
			rm $FIRST_PAIR_BAM

			printProgress "[projectAllelic stranded RPM] Splitting $FLAG_BAM reads by 2nd in pair..."
			local SECOND_PAIR_BAM=$USER_INPUT_NAME_MAP_FLAG"_secondInPair.bam"
			$SAMTOOLS view -bh -f 0x0080 $FLAG_BAM > $SECOND_PAIR_BAM
			$BEDTOOLS genomecov -ibam $SECOND_PAIR_BAM -bg -split -scale $SCALING_FACTOR -strand + > $USER_INPUT_NAME_MAP_FLAG"_second_pos_RPM.bedGraph"
			$BEDTOOLS genomecov -ibam $SECOND_PAIR_BAM -bg -split -scale $SCALING_FACTOR -strand - > $USER_INPUT_NAME_MAP_FLAG"_second_neg_RPM.bedGraph"
			rm $SECOND_PAIR_BAM

			#plus/pos strand
			printProgress "[projectAllelic stranded RPM] Combining stranded bedgraphs for plus strand..."
			$BEDTOOLS unionbedg -i $USER_INPUT_NAME_MAP_FLAG"_first_neg_RPM.bedGraph" $USER_INPUT_NAME_MAP_FLAG"_second_pos_RPM.bedGraph" > $USER_INPUT_NAME_MAP_FLAG"_p_tmp.bedGraph"
			rm $USER_INPUT_NAME_MAP_FLAG"_first_neg_RPM.bedGraph" $USER_INPUT_NAME_MAP_FLAG"_second_pos_RPM.bedGraph"
			awk '{OFS="\t";FS="\t"} {print $1, $2, $3, $4+$5}' $USER_INPUT_NAME_MAP_FLAG"_p_tmp.bedGraph" > $USER_INPUT_NAME_MAP_FLAG"_pos_preProject.bedGraph"
			rm $USER_INPUT_NAME_MAP_FLAG"_p_tmp.bedGraph"
			prepWigAndProject $USER_INPUT_NAME_MAP_FLAG"_RPM_pos" $USER_INPUT_NAME_MAP_FLAG"_pos_preProject.bedGraph" " plus stranded RPM"

			#minus/neg strand
			printProgress "[projectAllelic stranded RPM] Combining stranded bedgraphs for minus strand..."
			$BEDTOOLS unionbedg -i $USER_INPUT_NAME_MAP_FLAG"_first_pos_RPM.bedGraph" $USER_INPUT_NAME_MAP_FLAG"_second_neg_RPM.bedGraph" > $USER_INPUT_NAME_MAP_FLAG"_n_tmp.bedGraph"
			rm $USER_INPUT_NAME_MAP_FLAG"_first_pos_RPM.bedGraph" $USER_INPUT_NAME_MAP_FLAG"_second_neg_RPM.bedGraph"
			awk '{OFS="\t";FS="\t"} {print $1, $2, $3, $4+$5}' $USER_INPUT_NAME_MAP_FLAG"_n_tmp.bedGraph" > $USER_INPUT_NAME_MAP_FLAG"_neg_preProject.bedGraph"
			rm $USER_INPUT_NAME_MAP_FLAG"_n_tmp.bedGraph"
			prepWigAndProject $USER_INPUT_NAME_MAP_FLAG"_RPM_neg" $USER_INPUT_NAME_MAP_FLAG"_pos_preProject.bedGraph" " minus stranded RPM"

		else #unstranded RPM
			printProgress "[projectAllelic unstranded RPM] Converting $FLAG_BAM to bedGraph"
			$BEDTOOLS genomecov -ibam $FLAG_BAM -bg -split -scale $SCALING_FACTOR > $USER_INPUT_NAME_MAP_FLAG"_RPM_preProject.bedgraph"
			prepWigAndProject $USER_INPUT_NAME_MAP_FLAG"_RPM" $USER_INPUT_NAME_MAP_FLAG"_RPM_preProject.bedgraph" " unstranded RPM"
		fi
	done

	printProgress "[projectAllelic] All haplotype-specific bams have completed projection at [$(date)]"
	cd $CURRENT_DIRECTORY
}

function prepWigAndProject () {
	local FINAL_NAME=$1
	local PRE_BEDGRAPH=$2
	local PROGRESS_APPEND=$3

	printProgress "[projectAllelic$PROGRESS_APPEND] Converting $PRE_BEDGRAPH to WIG"
	local PRE_WIG=$FINAL_NAME"_preProject.wig"
	awk '
	BEGIN {
				 print "track type=wiggle_0"
	}
	NF == 4 {
					 print "fixedStep chrom="$1" start="$2+1" step=1 span=1"
					 for(i = 0; i < $3-$2; i++) {
									print $4
					}
	}' $PRE_BEDGRAPH > $PRE_WIG
	gzip $PRE_WIG

	printProgress "[projectAllelic$PROGRESS_APPEND ALEA] Projecting $PRE_WIG and refmap files to reference coordinates"
	local PROJECTED_BEDGRAPH=$FINAL_NAME".bedgraph" #output bedgraph with reference genome coordinates
	if [[ $PRE_WIG == *$HAPLO_1* ]]; then
		$ALEA project --input-wig=$PRE_WIG".gz" --input-refmaps=$HAPLO_1_REFMAP --output-bedgraph=$PROJECTED_BEDGRAPH
	elif [[ $PRE_WIG == *$HAPLO_2* ]]; then
		$ALEA project --input-wig=$PRE_WIG".gz" --input-refmaps=$HAPLO_2_REFMAP --output-bedgraph=$PROJECTED_BEDGRAPH
	fi

	mv $PROJECTED_BEDGRAPH $CURRENT_DIRECTORY/TrackHub/

	#maybe can integrate with masterTrackHub
	printProgress "[projectAllelic$PROGRESS_APPEND] Converting projected $PROJECTED_BEDGRAPH to BIGWIG"
	$BEDGRAPHTOBW $PROJECTED_BEDGRAPH $CHROM_SIZES $FINAL_NAME".bw"
}




# DEPRECATED BY JRA JULY 2020. RATIONALE:
# I think users should know what library type (ChIP, RNA, RRBS, PBAT, MethylSeq)
# each of the dataset is, and that they should include this info in the filename

# Drawbacks:
# - user friendliness decreased
# - input files or search keys must abide by a rigid rule (including library type)

# Positives: 
# - code decreased from 4 functions, 146 lines to 1 function, 54 lines
# - no longer rely on meta data provided by the authors, which is sometimes incorrect (including PAIRED-ENDEDNESS info)
# - sometimes uploading users input "OTHER" as library type

### Download reads from SRA
function downloadReads () {
	DL_COUNTER=0

	SRACODE=$code
	NAME=$(grep -e $SRACODE $INPUT_SAMPLE_LIST_FILE | cut -f2)

	printProgress "[masterDownload wget] Downloading $SRACODE to $TEMP_DIR... at [$(date)]"
	for DL in $($ESEARCH -db sra -query $SRACODE \
							| $EFETCH --format runinfo \
							| cut -d ',' -f 10 \
							| grep https); do
		wget --no-check-certificate $DL -O $SEARCH_KEY"_"$(basename $DL) & #& this allow the command to run in parallel and in the background 
		DL_PID_ARRAY[$DL_COUNTER]=$! #$! = the last process that was started
		((DL_COUNTER++)) #add one to the counter so next thing added to the array will be in the next position
	done

	for pid in ${DL_PID_ARRAY[*]}
	do
		wait $pid #wait for all downloads to finish
	done

	printProgress "[masterDownload] Finished $SRACODE -> $USER_INPUT_NAME reads download."
}



function extractType() {
	SRACODE=${1:-$SRACODE} #giving option to just use function to see what type of sequence it is? is this neccessary?

	TYPE=$($ESEARCH -db sra -query $SRACODE \
				| $EFETCH --format runinfo \
				| cut -d ',' -f 13 \
				| tail -n 2 \
				| head -n 1)
	printProgress "[masterDownload] Data type: $LIBRARY_TYPE"

	if [[ $TYPE == "ChIP-Seq" ]] && \
		 [[ $USER_INPUT_NAME != *"ChIP"* ]]; then
		 NAME="ChIPseq_"$USER_INPUT_NAME

	elif [[ $TYPE == "RNA-Seq" ]] && \
			 [[ $USER_INPUT_NAME != *"RNA"* ]] ; then
			 MIN_MAPQ=255
			 NAME="RNAseq_"$USER_INPUT_NAME

	elif [[ $TYPE == "Bisulfite-Seq" ]] && \
			 [[ $USER_INPUT_NAME != *"RRBS"* ]] && \
			 [[ $USER_INPUT_NAME != *"BSSeq"* ]] && \
			 [[ $USER_INPUT_NAME != *"PBAT"* ]] ; then
			 NAME="BSSeq_"$USER_INPUT_NAME

	elif [[ $TYPE == "ATAC-seq" ]] && \
			 [[ $USER_INPUT_NAME != *"ATAC"* ]]; then
			 NAME="ATACSeq_"$USER_INPUT_NAME

	elif [[ $TYPE == "DNase-Hypersensitivity" ]] && \
			 [[ $USER_INPUT_NAME != *"DNase"* ]]; then
			 NAME="DNase_"$USER_INPUT_NAME

	elif [[ $TYPE == "Hi-C" ]] && \
			 [[ $USER_INPUT_NAME != *"HiC"* ]]; then
			 NAME="HiC_"$USER_INPUT_NAME
	fi
	printProgress "[masterDownload] Name of data: $USER_INPUT_NAME"
}

### Determine data to be single or paired-end (could be called with just SRACODE)
function extractPaired () {
	SRACODE=${1:-$SRACODE}

	PAIRED=$($ESEARCH -db sra -query $SRACODE \
					| $EFETCH --format runinfo \
					| cut -d ',' -f 16 \
					| tail -n 2 \
					| head -n 1)

	if [[ $PAIRED == SINGLE ]]; then
		PAIRED_END=false
	printProgress "[masterDownload] Data are single-end."
	else
		PAIRED_END=true
	printProgress "[masterDownload] Data are paired-end."
	fi
}


### Extract fastq using fasterq-dump, compress FASTQ then concatenate FASTQ from same samples
function extractFastq () {
	COMPRESS_COUNTER=0
	COMPRESS_PID_ARRAY=""

	local SEARCH_DL=*$SEARCH_KEY*[DSE]RR*

	printProgress "[masterDownload extractFastq] Dumping fastq files..."
	for SRA_FILE in $SEARCH_DL; do
		$FASTERQDUMP -e $RUN_THREAD --split-files ./$SRA_FILE
	done

	printProgress "[masterDownload extractFastq] Compressing fastq files..." #simultaneously zip all the dump fastq files for that read
	for DL_FASTQ in $SEARCH_DL*fastq; do
		gzip $DL_FASTQ &
		COMPRESS_PID_ARRAY[$COMPRESS_COUNTER]=$!
		((COMPRESS_COUNTER++))
	done

	for gzipid in ${COMPRESS_PID_ARRAY[*]}; do
		wait $gzipid #wait for all gzip to finish
	done

	printProgress "[masterDownload extractFastq] Concatenating fastq.gz files..."
	if $PAIRED_END; then
		COUNT=$(ls -1 $SEARCH_DL*_1.fastq.gz | wc -l)
		if [[ $COUNT > 1 ]] ; then
		# JUlien added printProgress and checkFileExists calls here - untested so far (added 6 calls for each)
		printProgress "[masterDownload extractFastq] Concatenating $SEARCH_DL*_1.fastq.gz into $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME"_1.fastq.gz""
			cat $SEARCH_DL*_1.fastq.gz > $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME"_1.fastq.gz"
		printProgress "[masterDownload extractFastq] Concatenating $SEARCH_DL*_2.fastq.gz into $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME"_2.fastq.gz""
			cat $SEARCH_DL*_2.fastq.gz > $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME"_2.fastq.gz"
		checkFileExists $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME"_1.fastq.gz"
		checkFileExists $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME"_2.fastq.gz"
		else # If only single file, move, don't copy
		printProgress "[masterDownload extractFastq] Renaming $SEARCH_DL*_1.fastq.gz to $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME"_1.fastq.gz""
			mv $SEARCH_DL*_1.fastq.gz $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME"_1.fastq.gz"
		checkFileExists $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME"_1.fastq.gz"
		checkFileExists $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME"_2.fastq.gz"
		printProgress "[masterDownload extractFastq] Renaming $SEARCH_DL*_2.fastq.gz to $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME"_2.fastq.gz""
			mv $SEARCH_DL*_2.fastq.gz $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME"_2.fastq.gz"
		$CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME"_2.fastq.gz"
		fi

	else #Single-End Reads
		COUNT=$(ls -1 $SEARCH_DL*.fastq.gz | wc -l)
		if [[ $COUNT > 1 ]] ; then # If only single file, move, don't copy
			printProgress "[masterDownload extractFastq] Concatenating $SEARCH_DL*.fastq.gz into $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME".fastq.gz""
			cat $SEARCH_DL*.fastq.gz > $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME".fastq.gz"
			checkFileExists $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME".fastq.gz"
		else
			printProgress "[masterDownload extractFastq] Renaming $SEARCH_DL*.fastq.gz to $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME".fastq.gz""
			mv $SEARCH_DL*.fastq.gz $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME".fastq.gz"
			checkFileExists $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$USER_INPUT_NAME".fastq.gz"
		fi
	fi

	printProgress "[masterDownload] Fastq files for $USER_INPUT_NAME are moved to $CURRENT_DIRECTORY/$FASTQ_DRECTORY"
	rm $SEARCH_DL
}
