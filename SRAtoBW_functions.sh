#! /bin/bash

###List of functions that can be used in various pipelines for download, align and creating track hub

# TODO 2021-11-03 BAM files created under user group, not lab group on Graham
#DO THIS THING

##System Specific Configuration

#Only if you ever plan on using these functions directly in the command line (otherwise ignore)
#Ensure correct config file if using
#if [[ -z $SCRIPTS_DIR ]]; then
	# pushd $(dirname $0) > /dev/null
	# SCRIPTS_DIR=$(pwd -P)
	# popd > /dev/null
	# source $SCRIPTS_DIR/BRC.config
# fi


## Non System Specific Variables
CURRENT_DIRECTORY=$(pwd)
FASTQ_DIRECTORY="./Fastq"
TRACK_HUB_DIR="./Track_Hub"

ALLELE_SPECIFIC=false
ALLELE_RUN=false
BAM_INPUT=false
BAM_REALIGNMENT=false
FASTQ_INPUT=false
FASTQ_ONLY=false
KEEP_FASTQ=false
KEEP_REPLICATES=false
#LOCAL_RUN=false
#SEP_PARA=false
STRANDED_ALLELIC=false
TRIM_READ=false
USE_BOWTIE=false
CHECK_DEPEND=false

CODE_ARRAY=""
BIN_SIZE=1
FLAG=1540 #read unmapped, read fails platform/vendor quality checks, read is PCR or optical duplicate 
GENOME_BUILD="mm10"
MIN_MAPQ=5
NORMALIZE="CPM" # TODO Should this be BPM (aka TPM) for RNAseq data?
SMOOTH_WIN=0


# Help Menu
OPTIONS="hi:ab:B:d:Df:Fg:kK:Lm:M:n:N:ors:t:Tuw:x"

HELP="USAGE:\t $(basename $0) [OPTIONS] -h for help"

HELP_FULL="\n$HELP\n
\nThis script has been designed to streamline and standardize the procurement of data from the SRA database, as well as alignment and visualization of all data. Given an input file listing the SRA files to be downloaded, it will automatically download them, decompress, rename according to the users desires, align to the genome of choice, collapse replicates and finally generate a UCSC-compatible track hub including bigwig files. Users can opt to use some or all of these features depending on their input or output needs.\n
\nIt is generally best to run this script in a clean folder, as it performs many cleanup steps on both the current folder and ALL subdirectories. For instance, any folders that include \"STAR\" in the name will be removed, as the STAR program leaves temporary aligns behind when finished.\n
\nAlignment:\tChIPseq - bwa mem\n\t\tRNAseq - STAR\n\t\tBSSeq/PBAT/RRBS - bismark\n
\n\nOPTIONS:
-h\tPrints help page.\n\t
-i\tSRA Input File. Must be in the format Accession(tab)DesiredName.\n\t\tDesiredName must be formatted specifically. A final identifier \n\t\tmust be present at the end of the name by which data should \n\t\tbe grouped (ex. Foo2018, etc.). The aligned files will be in a \n\t\tdirectory of this name. If there are replicates, the name will\n\t\tbe followed by an underscore and \"RepX\" - this will be removed \n\t\tfrom the final name upon collapsing replicates. Data type \n\t\tshould be included somewhere in the name; if not, it will be\n\t\tprepended to the name as 'RNAseq' or 'BSSeq' or 'ChIPseq'.\n\t
-a\tAllele-specific alignment using MEA.\n\t
-b\tBAM Inputs. Will use any already aligned .bam files in\n\t\tsubdirectories to generate trackhub.\n\t
-B\tBAM Input with realignment. Will extract reads from .bam files\n\t\tin listed directory and realign to genome of choice.\n\t
-d\tTemporary directory. Useful for solid-state drives etc.\n\t\tPlease change the default using BRC.config.\n\t
-D\tCheck Dependencies and exit.\n\t
-f\tInput .fastq files. Provide folder name in which fastq files\n\t\tare located (files must end in .fastq.gz).\n\t
-F\tOnly output .fastq files.\n\t
-g\tGenome build for alignment. Allowable: mm9, mm10, rn5, rn6,\n\t\toryCun2, mesAur1, or hg19. Default=mm10\n\t
-k\tKeep .fastq files when done.\n\t
-K\tNCBI API Key for use.\n\t
-m\tMemory to give each thread (Format=XG/M).\n\t
-M\tMinimum mapping quality for bigwig generation. Default=5\n\t
-n\tBin size for bigwig generation. Larger bins to smooth noisy\n\t\tdata. Default=1\n\t
-N\tNormalization method for bigwigs. Accepted: CPM, RPKM\n\t\t(Default=CPM)\n\t
-r\tKeep replicates after collapsing. Default=false.\n\t
-s\tObtained stranded RPM tracks for allele-specific runs.\n\t\tDefault=false\n\t
-t\tNumber of Threads to use. Default=6 (Check in config file)\n\t
-T\tTrim .fastq files after download.\n\t
-u\tChanging ChIPseq aligner to bowtie2. Default=BWA\n\t
-w\tSmoothing window. Will smooth bigwigs in a rolling window of\n\t\tthis size. Default=0\n\t"

#-L\tRunning script on a local computer rather than a server.\n\t
#-x\tSpecify search key of data to use if want to apply pipeline on \n\t\tspecific set of data when entering pipeline through fastq or bam.\n\t
#-X\t?\n\t

#TODO -x option; -X option?!?



########################################
#							 FUNCTIONS							 #
########################################

### Changing options/variables based on options passed from shell script
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
#				SEP_PARA=true
#				PASS_ARG=$@
				FILENAME=${OPTARG}
				;;
			a)
				ALLELE_SPECIFIC=true
				;;
			b) #bam input for trackhub
				BAM_INPUT=true
				KEEP_REPLICATES=true
				if $FASTQ_ONLY ; then
					echo -e "ERROR:\tIncompatible options - track hub generation downstream of bam input." 
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
				CHECK_DEPEND=true
				;;
			f) #use existing fastq files for downstream modules in pipeline (will skip download)
				FASTQ_INPUT=true
				KEEP_FASTQ=true
				if $FASTQ_ONLY ; then
					echo -e "ERROR:\tFiles already in Fastq format."
					exit 1
				fi
				FASTQ_DIRECTORY=${OPTARG}
				;;
			F) #Stop pipeline after downloading
				FASTQ_ONLY=true
				if $BAM_INPUT ; then
					echo -e "ERROR:\tIncompatible options. Cannot generate .fastq files from .bam files." # TODO this definitely isn't true
					exit 1
				fi
				;;
			g)
				GENOME_BUILD=${OPTARG}
				;;
			k) #keep fastq files after alignment
				KEEP_FASTQ=true
				;;
			K)
				export NCBI_API_KEY=${OPTARG}
				;;
#			L)
#				LOCAL_RUN=true
#				;;
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
#				if [[ $RUN_THREAD > 1 ]] ; then
#					((RUN_THREAD--))
#				fi
				;;
			T) #trimming fastq by default parameters
				TRIM_READ=true
				;;
			u)
				USE_BOWTIE=true
				;;
			w)
				SMOOTH_WIN=${OPTARG}
				;;
#			x)
#				SEP_PARA=false
#				SEARCH_KEY=${OPTARG}
#				;;
#			X)
#				CODE_ARRAY=$(echo $@ | sed 's/.*-X //') #TODO wtf is this doing ANSWER: apparently repopulating CODE_ARRAY give an input array.
#				;;
			\?)
				echo -e "\n###############\nERROR: Invalid Option! \nTry '$(basename $0) -h' for help.\n###############" >&2
				exit 1
				;;
		esac
	done

#	if $SEP_PARA; then
#		checkDependencies
#		setGenome $GENOME_BUILD
#	fi
}

### Checking dependencies of the functions TODO not sure this works with java dependencies properly
function checkDependencies () {
	DEPENDENCIES=("$TRIMMOMATIC" $STAR $BISMARK $BOWTIE2 $BWA $SAMTOOLS "$PICARD" awk $BAM2FASTQ $BEDGRAPHTOBW $BAMCOVERAGE)

	echo -e "[checkDependencies] Checking Dependencies"
	EXIT=0
	for COMMAND in "${DEPENDENCIES[@]}"; do
		echo -e "[checkDependencies] $COMMAND..."
		command -v $COMMAND > /dev/null 2>&1 || {
			echo -e >&2 "\t\t$COMMAND not found!"
			EXIT=1
		}
	done

	if [[ $EXIT = 1 ]] ; then
		exit 1
	fi
}

###Set up log file for runs
function setUp () {
#	if [[ $SEP_PARA == false ]]; then 
	#	if [[ -z $SEARCH_KEY ]]; then
			#without a search term associated with a -f or -b command, the logfiles might be messy (written to same logfile)
			#using a random number in placement to avoid confusion
	#		LOG_FILE=$CURRENT_DIRECTORY/$RANDOM"_"$(date '+%y-%m-%d')"_log.txt"
	#	else
			LOG_FILE=$CURRENT_DIRECTORY/$(date '+%y-%m-%d')"_log.txt"
	#	fi
	

		checkDependencies		
		if [[ $CHECK_DEPEND == "true" ]]; then
			exit
		fi

		printProgress "[setUp] Starting Script"
#		printProgress "[setUp] Search key for the set: $SEARCH_KEY"
		printProgress "[setUp] SRA array: $CODE_ARRAY"
		printProgress "[setUp] All required dependencies are found"
		setGenome $GENOME_BUILD
#	fi
}

function printProgress () {
	echo -e $(date '+%Y-%b-%d %T') $1 | tee -a $LOG_FILE #using tee will also show in stdout
}

function checkFileExists () {
	if [[ ! -f $1 ]]; then
		echo -e "ERROR:\tFile $1 does not exist"
		exit 1
	fi
}

### Create neccessary files for reference genome
function setGenome () {
#	GENOME_BUILD=$1

	#check that the directory exist here
	if [[ -d $GENOME_DIR/$GENOME_BUILD ]]; then
			GENOME_DIR=$GENOME_DIR/$GENOME_BUILD
#			if [[ $SEP_PARA == false ]]; then
				printProgress "[setGenome] Genome used for data is $GENOME_BUILD"
				printProgress "[setGenome] Path to genome directory: $GENOME_DIR"
#			fi
	else
		echo -e "Genome entered is invalid, please check reference genome folder."
		exit
	fi


#Specify location of the following files & directories if layout is not as example
#The $GENOME_DIR specified in this location, already includes the specific species
	CHROM_SIZES=$GENOME_DIR/$GENOME_BUILD".sizes"
	GENOME_FILE=$GENOME_DIR/$GENOME_BUILD".fa"
	ACTB_BED=$GENOME_DIR/$GENOME_BUILD"_Actb.bed"

	STAR_GENOME_DIR=$GENOME_DIR/$GENOME_BUILD"-STAR"
	STAR_SJ_DB=$GENOME_DIR/$GENOME_BUILD"-NCBIRefSeq.gtf"
	BISMARK_GENOME_DIR=$GENOME_DIR # TODO should this be /bisulfite?
	BOWTIE2_INDEXES=$GENOME_DIR/$GENOME_BUILD

	#directory that stores pseudogenomes for allelic pipeline
	DIPLOID_GENOME_DIR=$GENOME_DIR/diploid
	HAPLOID_GENOME_DIR=$GENOME_DIR/haploid
	
	checkFileExists $CHROM_SIZES
	checkFileExists $GENOME_FILE
	checkFileExists $ACTB_BED




}

### Create subset of SRACODE array to be used in parallel running
# TODO Why is this so complicated and opaque?
# function parallelRun () {
	# if $SEP_PARA; then
		# echo -e "Separating input file into subsets for parallel runs..."

		# INPUT_FILE=$CURRENT_DIRECTORY/$FILENAME
		# createCodeArray
		# CURRENT_SET=""
	
		# for code in $CODE_ARRAY; do
			# SRACODE=$code
			# NAME=$(grep -e $SRACODE $INPUT_FILE | cut -f2)

			# if [[ $NAME == *_[Rr]ep* ]] ; then
				# if [[ $CURRENT_SET != ${NAME//_[Rr]ep*/}* ]]; then
					# CURRENT_SET=${NAME//_[Rr]ep*/}
					# declare -a SUB_ARRAY=$(grep -e $CURRENT_SET $INPUT_FILE | cut -f1)
					# echo "calling "$(basename $SHELL_SCRIPT) "on" $CURRENT_SET
					
					# if $LOCAL_RUN; then
						# $SHELL_SCRIPT $PASS_ARG -x $CURRENT_SET -X $SUB_ARRAY & 
						# wait $! #wait for the script above to finish running before moving onto the next set (avoid overload)
					# else
						# echo "Submitting on server"
						# DATE=$(date '+%y-%m-%d')
						# $SERVER_SUBMIT "MasterDAT_"$DATE"_"$CURRENT_SET $SHELL_SCRIPT $PASS_ARG -x $CURRENT_SET -X $SUB_ARRAY
						# sleep 10 #pause before running next code 
					# fi

				# fi

			# else
				# declare -a SUB_ARRAY=$SRACODE
				# echo "calling "$(basename $SHELL_SCRIPT) "on" $NAME

				# if $LOCAL_RUN ; then
					# $SHELL_SCRIPT $PASS_ARG -x $NAME -X $SUB_ARRAY
					# wait $!
				# else
					# echo "Submitting on server"
					# DATE=$(date '+%y-%m-%d')
					# $SERVER_SUBMIT "MasterDAT_"$DATE"_"$NAME $SHELL_SCRIPT $PASS_ARG -x $NAME -X $SUB_ARRAY
					# sleep 10
				# fi
#
#			fi
#		done
#		exit #exit the script b/c don't want to run the rest of the code on every single SRACODE again
#	fi

#}

### Downloading files specified from the tab-delimited file to fastq files
function masterDownload () {
	if $BAM_INPUT || $FASTQ_INPUT; then #if input is either BAM/FASTQ then no download is needed - exit masterDownload function
		return 0
	fi

	mkdir $FASTQ_DIRECTORY
	if $BAM_REALIGNMENT; then
		extractFastqFromBAM
		return 0 #after extracting fastq from BAM, can exit masterDownload function
	fi

	INPUT_FILE=$CURRENT_DIRECTORY/$FILENAME
	
### Create an array of SRACODE passed listed in the tab-delimited file

	while read n; do #read command will read a line of and split them into words (2 words in files)
		declare -a m=($n) #-a for an array m to be stored based on the line (SRACODE, NAME)
		CODE_ARRAY=$CODE_ARRAY" "${m[0]} #add only SRACODE to this list
	done < $INPUT_FILE #feeding FILE as stdin to this while read

	
	cd $TEMP_DIR

	printProgress "[masterDownload] Starting..."

	for CODE in $CODE_ARRAY; do
		downloadReads $CODE
		$ESEARCH -db sra -query $CODE | $EFETCH -format runinfo > "temp.txt" # Reduce the number of times we query NCBI servers
		extractType
		extractPaired
		extractFastq $CODE
	done
	rm "temp.txt"
	cd $CURRENT_DIRECTORY

	printProgress "[masterDownload] All fastq files are downloaded."

	if $FASTQ_ONLY; then
		printProgress "[masterDownload] Fastq files only. Exit script."
		exit 0 #exit the whole script b/c only requires fastq files
	fi
}
 



### Download reads from SRA
function downloadReads () {
	DL_COUNTER=0

	SRACODE=$1
	NAME=$(grep -e $SRACODE $INPUT_FILE | cut -f2)

#	sleep $(( $RANDOM % 300 )) #sleep for random (between 0 and 300) seconds TODO WHYYyyyyy
	printProgress "[masterDownload wget] Downloading $SRACODE to $TEMP_DIR..."
#	for DL in $($ESEARCH -db sra -query $SRACODE \ 
#							| $EFETCH -format runinfo \
#							| cut -d ',' -f 10 \
#							| grep https); do			
#		wget -q --no-check-certificate $DL & #& this allow the command to run in parallel and in the background  TODO might have to change this to prefetch
	$PREFETCH -X 100G -O $TEMP_DIR $SRACODE &
	DL_PID_ARRAY[$DL_COUNTER]=$! #$! = the last process that was started
	((DL_COUNTER++)) #add one to the counter so next thing added to the array will be in the next position
#	done
		
	for pid in ${DL_PID_ARRAY[*]}
	do
		wait $pid #wait for all downloads to finish
	done

	printProgress "[masterDownload] Finished $SRACODE -> $NAME reads download."
}

### Determine data type (could be called with just SRACODE)
function extractType() {
	SRACODE=${1:-$SRACODE} #giving option to just use function to see what type of sequence it is? is this neccessary?

	TYPE=$(cat "temp.txt" \
				| cut -d ',' -f 13 \
				| head -n 2 \
				| tail -n 1)
	printProgress "[masterDownload] Data type: $TYPE"

	if [[ $TYPE == "ChIP-Seq" ]] && \
		 [[ $NAME != *"ChIP"* ]]; then
		NAME="ChIPseq_"$NAME

	elif [[ $TYPE == "RNA-Seq" ]] && \
			 [[ $NAME != *"RNA"* ]] ; then
			 MIN_MAPQ=255
			 NAME="RNAseq_"$NAME

	elif [[ $TYPE == "Bisulfite-Seq" ]] && \
			 [[ $NAME != *"RRBS"* ]] && \
			 [[ $NAME != *"BSSeq"* ]] && \
			 [[ $NAME != *"PBAT"* ]] ; then
			 NAME="BSSeq_"$NAME

	elif [[ $TYPE == "ATAC-seq" ]] && \
			 [[ $NAME != *"ATAC"* ]]; then
			 NAME="ATACSeq_"$NAME

	elif [[ $TYPE == "DNase-Hypersensitivity" ]] && \
			 [[ $NAME != *"DNase"* ]]; then
			 NAME="DNase_"$NAME

	elif [[ $TYPE == "Hi-C" ]] && \
			 [[ $NAME != *"HiC"* ]]; then
			 NAME="HiC_"$NAME
	fi
	printProgress "[masterDownload] Name of data: $NAME"
}

### Determine data to be single or paired-end (could be called with just SRACODE)
function extractPaired () {
	SRACODE=${1:-$SRACODE}

	PAIRED=$(cat "temp.txt" \
					| cut -d ',' -f 16 \
					| head -n 2 \
					| tail -n 1)

	if [[ $PAIRED == SINGLE ]]; then
		PAIRED_END=false
	printProgress "[masterDownload] Data are single-end."
	else
		PAIRED_END=true
	printProgress "[masterDownload] Data are paired-end."
	fi
} 

### Extract fastq using fasterq-dump, compressed then concatenate fastq
function extractFastq () { 
	COMPRESS_COUNTER=0
	COMPRESS_PID_ARRAY=""
	#CODE=$1

#	local SEARCH_DL="*[DSE]RR*"

	printProgress "[masterDownload extractFastq] Dumping fastq files..."
	for SRA_FILE in */*.sra; do
		$FASTERQDUMP -e $RUN_THREAD --split-files $SRA_FILE
		rm $SRA_FILE
	done

	printProgress "[masterDownload extractFastq] Compressing fastq files..." #simultaneously zip all the dump fastq files for that read
	for DL_FASTQ in *.fastq; do 
		gzip $DL_FASTQ & 
		COMPRESS_PID_ARRAY[$COMPRESS_COUNTER]=$!
		((COMPRESS_COUNTER++))
	done

	for gzipid in ${COMPRESS_PID_ARRAY[*]}; do
		wait $gzipid #wait for all gzip to finish
	done

	printProgress "[masterDownload extractFastq] Concatenating fastq.gz files..."
	if $PAIRED_END; then
		COUNT=$(ls -1 *_1.fastq.gz | wc -l)
		if [[ $COUNT > 1 ]] ; then 
			cat *_1.fastq.gz > $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_1.fastq.gz"
			cat *_2.fastq.gz > $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_2.fastq.gz"
		else # If only single file, move, don't copy
			mv *_1.fastq.gz $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_1.fastq.gz"
			mv *_2.fastq.gz $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_2.fastq.gz"
		fi

		#checkpoint: fastq files have been moved correctly
		checkFileExists $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_1.fastq.gz"
		checkFileExists $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_2.fastq.gz"
 
	else #Single-End Reads
		COUNT=$(ls -1 *.fastq.gz | wc -l)
		if [[ $COUNT > 1 ]] ; then # If only single file, move, don't copy
			cat *.fastq.gz > $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME".fastq.gz"
		else
			mv *.fastq.gz $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME".fastq.gz"
		fi

		#checkpoint: fastq files have been moved correctly
		checkFileExists $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME".fastq.gz"

	fi

	printProgress "[masterDownload] Fastq files for $NAME are moved to $CURRENT_DIRECTORY/$FASTQ_DIRECTORY"
	rm *.fastq.gz
	rm -r $CODE
}

### Trim fastq files for adaptors using trimmomatic
### Could be called with $1=directory that holds fastq files TODO Could it?
function trimReads () {
	if [[ $TRIM_READ == false || -z $1 ]]; then #exit script if not needed
		printProgress "[trimReads] Reads not trimmed."
		return 0
	fi

	cd $TEMP_DIR
	printProgress "[trimReads] Starting..."
		
	for FILE in $CURRENT_DIRECTORY/${1:-$FASTQ_DIRECTORY}/*fastq.gz; do
		determinePairedFastq

		if [[ $PAIRED_END ]]; then
			printProgress "[trimReads] Trimming "$NAME"_1.fastq.gz and "$NAME"_2.fastq.gz..."

			$TRIMMOMATIC PE -threads 10 $FASTQ_PATH"_1.fastq.gz" $FASTQ_PATH"_2.fastq.gz" \
			$NAME"_trim_1.fastq.gz" $NAME"_unpaired_trim_1.fastq.gz" $NAME"_trim_2.fastq.gz" $NAME"_unpaired_trim_2.fastq.gz" \
			ILLUMINACLIP:$ILLUMINA_ADAPATORS_ALL":2:30:10" \
			SLIDINGWINDOW:4:20 \
			MINLEN:36

			#checkpoint: fastq have been trimmed correctly
			checkFileExists $NAME"_trim_1.fastq.gz"
			checkFileExists $NAME"_trim_2.fastq.gz"
			
			mv $NAME"_trim_1.fastq.gz" $FASTQ_PATH"_1.fastq.gz"
			mv $NAME"_trim_2.fastq.gz" $FASTQ_PATH"_2.fastq.gz"

			#checkpoint: fastq files have been moved correctly
			checkFileExists $FASTQ_PATH"_1.fastq.gz"
			checkFileExists $FASTQ_PATH"_2.fastq.gz"

		else #Single-End
			printProgress "[trimReads] Trimming $NAME.fastq.gz"

			$TRIMMOMATIC SE -threads 10 $FASTQ_PATH".fastq.gz" $FASTQ_PATH"_trim.fastq.gz" \
			ILLUMINACLIP:$ILLUMINA_ADAPATORS_ALL":2:30:10" \
			SLIDINGWINDOW:4:20 \
			MINLEN:36

			#checkpoint: fastq have been trimmed correctly
			checkFileExists $NAME"_trim.fastq.gz"

			mv $NAME"_trim.fastq.gz" $FASTQ_PATH".fastq.gz"

			#checkpoint: fastq files have been moved correctly
			checkFileExists $FASTQ_PATH".fastq.gz"
		fi

		printProgress "[trimReads] Finished trimming and renaming fastq files for $NAME..."
			
	done
		
	rm *unpaired*
	cd $CURRENT_DIRECTORY
}
	
### Determine paired or single-end from name of fastq file
function determinePairedFastq () {
	if [[ $FILE == *"_2.fastq.gz" ]] ; then 
		continue #dont process 2nd read, continue to next iteration
	fi
	
	if [[ $FILE == *"_1.fastq.gz" ]] ; then
		FASTQ_PATH=${FILE//_1.fastq.gz/} #removes everything that's after // - leaving path to directory
		NAME=${FASTQ_PATH##*/} #removes all prefixes prior to the last / - removing path to directory
		
		PAIRED_END=true
		FILE_FASTQ1=$FASTQ_PATH"_1.fastq.gz"
		FILE_FASTQ2=$FASTQ_PATH"_2.fastq.gz"

		#checkpoint: fastq files are found before passing into aligner
		checkFileExists $FILE_FASTQ1
		checkFileExists $FILE_FASTQ2

	else
		FASTQ_PATH=${FILE//.fastq.gz}
		NAME=${FASTQ_PATH##*/}

		PAIRED_END=false
		FILE_FASTQ=$FASTQ_PATH".fastq.gz"

		#checkpoint: fastq files are found before passing into aligner
		checkFileExists $FILE_FASTQ
	fi

	FILE_RAW_BAM=$NAME"_raw.bam"
	FILE_BAM=$NAME".bam"
	
	if [[ $NAME == *_[Rr]ep* ]] ; then #creating name of folder that will be storing the BAMs
		X=${NAME%_*} #removing the "_Rep"
		BAM_FOLDER_NAME=${X##*_} #Removing everything before the last _ (leaving grouping identifier)
	else
		BAM_FOLDER_NAME=${NAME##*_}
	fi

	BAM_FOLDER=$CURRENT_DIRECTORY/$BAM_FOLDER_NAME
	mkdir -p $BAM_FOLDER
}

### Align fastq files to using data-specific aligner
function masterAlign () {
	if $BAM_INPUT; then #if input is BAM for trackhub - exit function
		return 0
	fi

	cd $TEMP_DIR
	printProgress "[masterAlign] Starting..."
	
	STAR_ARGUMENTS="--genomeDir $STAR_GENOME_DIR --sjdbGTFfile $STAR_SJ_DB --runThreadN $RUN_THREAD --sjdbOverhang 70 --outFilterType BySJout --twopassMode Basic --twopass1readsN 1000000000 --outSAMunmapped Within --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --readFilesCommand zcat "

	#here kinda assumes that samtools and bowtie2 are both in the path, or else you will need to specify their path here
	BISMARK_ARGUMENTS="--temp_dir $TEMP_DIR --gzip --local --parallel $BISMARK_THREAD --bowtie2 --bam $BISMARK_GENOME_DIR"
#	SEARCH_KEY=${1:-$SEARCH_KEY}

	for FILE in $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/*fastq.gz; do
		determinePairedFastq # Assigns Fastq names

		if [[ $FILE == *"RNA"* ]]; then
			alignSTAR
		elif [[ $FILE == *"RRBS"* ]] || [[ $FILE == *"BSSeq"* ]] || [[ $FILE == *"PBAT"* ]]; then
			alignBismark
		elif [[ $FILE == *"HiC" ]]; then
			alignHiCUP
			continue #for HiC data, you don't need to refine BAM
		elif $ALLELE_SPECIFIC || $USE_BOWTIE; then
			alignBowtie2
		else
			alignBWA
		fi

		checkFileExists $FILE_RAW_BAM

		if $ALLELE_RUN; then #allele specific run will need to do unpacking
			continue #move onto next set instead of refining the produced RAW BAM
		fi
		
		refineBam
		
	done

	printProgress "[masterAlign] Alignment of fastq files to $GENOME_BUILD completed."

	cd $CURRENT_DIRECTORY

	if [[ $TYPE == "Hi-C" ]]; then
		exit #for HiC data, the rest of script doesn't need to be run (only to align with download and align with HiCUP"\
	fi
	
}

### Alignment of RNASeq data using STAR
function alignSTAR () {
	FILE_STAR_OUTPUT=$NAME"Aligned.out.bam"

	if $PAIRED_END; then
		printProgress "[masterAlign STAR] Aligning "$NAME"_1.fastq.gz and "$NAME"_2.fastq.gz to genome..."
		$STAR $STAR_ARGUMENTS --readFilesIn $FILE_FASTQ1 $FILE_FASTQ2 --outFileNamePrefix $NAME

	else #Single-End
		printProgress "[masterAlign STAR] Aligning $NAME.fastq.gz to genome..."
		$STAR $STAR_ARGUMENTS --readFilesIn $FILE_FASTQ --outFileNamePrefix $NAME
	fi

	checkFileExists $FILE_STAR_OUTPUT

	printProgress "[masterAlign STAR] Alignment completed -> $FILE_RAW_BAM"
	mv $FILE_STAR_OUTPUT $FILE_RAW_BAM

	mv $NAME"Log.final.out" $CURRENT_DIRECTORY
	rm $NAME*"out"*
	rm -r *"STAR"*
}

### Alignment of BSSeq data using Bismark
function alignBismark () {
	if [[ $NAME == *"PBAT"* ]] ; then
		BISMARK_ARGUMENTS="--pbat "$BISMARK_ARGUMENTS
	else
		BISMARK_ARGUMENTS="--non_directional "$BISMARK_ARGUMENTS
	fi

	if $PAIRED_END; then
		printProgress "[masterAlign Bismark] Aligning "$NAME"_1.fastq.gz and "$NAME"_2.fastq.gz to genome..."
		$BISMARK $BISMARK_ARGUMENTS -1 $FILE_FASTQ1 -2 $FILE_FASTQ2
		BISMARK_OUTPUT=$(basename $FILE_FASTQ1)"_bismark_bt2_pe.bam"

	else #Single-End
		printProgress "[masterAlign Bismark] Aligning $NAME.fastq.gz to genome..."
		$BISMARK $BISMARK_ARGUMENTS $FILE_FASTQ
		BISMARK_OUTPUT=$(basename $FILE_FASTQ)"_bismark_bt2.bam"
	fi

	printProgress "[masterAlign Bismark] Alignment completed -> $FILE_RAW_BAM"
	BISMARK_OUTPUT=${BISMARK_OUTPUT//.fastq.gz/}
	mv $BISMARK_OUTPUT $FILE_RAW_BAM
	rm *report.txt
}

### Alignment of HiC data with HiCUP, assumed to be paired end 
function alignHiCUP () {
	ENZYME=$(grep -e $NAME $INPUT_FILE | cut -f5 | head -n 1) 

	if [[ -z $ENZYME ]]; then
		printProgress "[masterAlign HiCUP] ERROR:\tDigestion enzyme was not entered for HiC dataset..."
		exit 1
	fi
		
	printProgress "[masterAlign HiCUP] Digesting $GENOME_BUILD with $ENZYME"
	hicup_digester --re1 $ENYZME --outdir $TEMP_DIR --genome $GENOME_BUILD $GENOME_FILE
	local DIGEST_FILE=$(ls $TEMPDIR/Digest*)

	printProgress "[masterAlign HiCUP] Aligning $NAME to $GENOME_BUILD"
	$HICUP --zip --bowtie2 $BOWTIE2 --digest $DIGEST_FILE --index $GENOME_DIR/$GENOME_BUILD --outdir $BAM_FOLDER --temp $TEMP_DIR --threads $RUN_THREAD $NAME"_1.fastq.gz" $NAME"_2.fastq.gz"

	printProgress "[masterAlign HiCUP] Alignment completed -> $FILE_RAW_BAM"
	rm $DIGEST_FILE
}

### Alignment of ChIP data with BWA
function alignBWA () {
	FILE_SAM=$NAME".sam"
		 
	if $PAIRED_END; then
		printProgress "[masterAlign BWA] Aligning "$NAME"_1.fastq.gz and "$NAME"_2.fastq.gz to genome..."
		$BWA mem -t $RUN_THREAD $GENOME_FILE $FILE_FASTQ1 $FILE_FASTQ2 > $FILE_SAM
	
	else #Single-End
		printProgress "[masterAlign BWA] Aligning "$NAME" to genome..."
		$BWA mem -t $RUN_THREAD $GENOME_FILE $FILE_FASTQ > $FILE_SAM 
	fi

	printProgress "[masterAlign BWA] Converting sam file to .bam file..."
	$SAMTOOLS view -bhS -@ $RUN_THREAD $FILE_SAM > $FILE_RAW_BAM
	printProgress "[masterAlign BWA] Alignment completed -> $FILE_RAW_BAM"
	rm $FILE_SAM
	
}

### Alignment of ChIPSeq data using BT2 (for allelic specific)
function alignBowtie2 () {
	FILE_SAM=$NAME".sam"

	if $PAIRED_END; then
		printProgress "[masterAlign Bowtie2] Aligning "$NAME"_1.fastq.gz and "$NAME"_2.fastq.gz to genome..."
		$BOWTIE2 -x $BOWTIE2_INDEXES --local -p $RUN_THREAD -1 $FILE_FASTQ1 -2 $FILE_FASTQ2 -S $FILE_SAM

	else #Single-End
		printProgress "[masterAlign Bowtie2] Aligning $NAME.fastq to genome..."
		$BOWTIE2 -x $BOWTIE2_INDEXES --local -p $RUN_THREAD -U $FILE_FASTQ -S $FILE_SAM
	fi

	printProgress "[masterAlign Bowtie2] Converting sam file to .bam file..."
	$SAMTOOLS view -bhS -@ $RUN_THREAD $FILE_SAM > $FILE_RAW_BAM
	printProgress "[masterAlign Bowtie2] Alignment completed -> $FILE_RAW_BAM"
	rm $FILE_SAM

}

### Sort and mark duplicates in BAM file after alignment
function refineBam () {
	local FILE_CLEANED_BAM=$NAME"_cleaned.bam"
	local FILE_SORTED_BAM=$NAME"_sort.bam"
	
	#soft clipping alignment that hangs off end of reference & set MAPQ to 0 for unmapped reads
	printProgress "[refineBAM] Refining $FILE_RAW_BAM"
	$PICARD CleanSam I=$FILE_RAW_BAM O=$FILE_CLEANED_BAM

	printProgress "[refineBAM] Sorting by coordinates..."
	$SAMTOOLS sort -@ $RUN_THREAD -m $SORT_MEM -o $FILE_SORTED_BAM -T $NAME $FILE_CLEANED_BAM
		
	printProgress "[refineBAM] Marking duplicates..." #not removing the duplicates
	$PICARD MarkDuplicates I=$FILE_SORTED_BAM O=$FILE_BAM M=$NAME"_markDupeMetrics.txt" TMP_DIR=$TEMP_DIR

	chgrp $GROUP $FILE_BAM
	mv $FILE_BAM $BAM_FOLDER
	checkFileExists $BAM_FOLDER/$FILE_BAM
	printProgress "[refineBAM] Final $FILE_BAM is moved to $BAM_FOLDER."
	
	rm $FILE_RAW_BAM $FILE_CLEANED_BAM $FILE_SORTED_BAM #remove all the buffer bam files
	rm *.txt

}

### Searches for .bam files containing _Rep#.bam in nested directories and combines them into a single bam #TODO make this more specific
function collapseReplicates () {
	cd $TEMP_DIR
	CURRENT=""
	printProgress "[collapseReplicates] Starting..."
	
	for FILE in $CURRENT_DIRECTORY/*/*.bam; do
		if [[ $FILE == *_[Rr]ep* ]]; then 
			MERGED_BAM=${FILE//_[Rr]ep*.bam/.bam}
			
			if [[ $CURRENT != ${FILE//_[Rr]ep*.bam/}*.bam ]]; then 
				local REP_COUNT=$(ls -l ${FILE//_[Rr]ep*.bam/}*.bam | wc -l)
				CURRENT=$FILE
				printProgress "[collapseReplicates] Merging $REP_COUNT replicates to "$(basename $MERGED_BAM)" and indexing"
				$SAMTOOLS merge --threads $RUN_THREAD --write-index -f $MERGED_BAM ${FILE//_[Rr]ep*.bam/}*.bam
				
				if [[ $KEEP_REPLICATES == false ]]; then
					printProgress "[collapseReplicates] Removing replicates"
					rm ${FILE//_[Rr]ep*.bam/}_*ep*.bam
				fi

#				printProgress "[collapseReplicates] Indexing BAM file..."
#				$SAMTOOLS index ${MERGED_BAM//.bam/}*.bam #index merged & replicates (if available)

				if $KEEP_REPLICATES; then
					local REP_DIR=$(dirname $FILE)/"Reps" 
					mkdir -p $REP_DIR
					printProgress "[collapseReplicates] Moving all replicates into $REP_DIR"
					mv ${MERGED_BAM//.bam/}*_[Rr]ep* $REP_DIR 
				fi
			fi
			
		else
			printProgress "[collapseReplicates] No replicates for $FILE"
			MERGED_BAM=$FILE
			printProgress "[collapseReplicates] Indexing BAM file..."
			$SAMTOOLS index $MERGED_BAM
		fi

		checkFileExists $MERGED_BAM

		printProgress "[collapseReplicates] Obtaining flagstats for "$(basename $MERGED_BAM)" flagstats"
		$SAMTOOLS flagstat $MERGED_BAM | tee -a $LOG_FILE
		
	done

	printProgress "[collapseReplicates] All replicates have been combined."

	cd $CURRENT_DIRECTORY
}

###Remove downloaded fastq files unless specified
# TODO what if this doesn't exist?
function removeFASTQ () {
	if [[ $KEEP_FASTQ == false ]]; then
		echo "Not keeping fastq..."
		local FASTQ_DIR=$CURRENT_DIRECTORY/$FASTQ_DIRECTORY
		rm -r $FASTQ_DIR
	fi
}

###Extract fastq file(s) from all bam within provided directory
function extractFastqFromBAM () {
	cd $TEMP_DIR
	printProgress "[extractFastqFromBAM] Starting..."
	
	for BAM_INPUT in $CURRENT_DIRECTORY/${1:-$BAM_REALIGNMENT_DIRECTORY}/*.bam; do
		BAM_NAME=$(basename $BAM_INPUT)
		local TEMP="temp"
		local SORTED="temp_sorted_"$BAM_NAME
		mkdir -p $TEMP
		OUTPUT=$TEMP/${BAM_NAME//.bam/#.fastq} #the # will be replaced by _1/_2 for PE reads in bam2fastq

		printProgress "[extractFastqFromBAM] Sorting $BAM_NAME by read name..."
		$SAMTOOLS sort -@ $RUN_THREAD -m $SORT_MEM -o $SORTED -n $BAM_INPUT

		printProgress "[extractFastqFromBAM] Extracting reads from $SORTED..."
		$BAM2FASTQ -o $OUTPUT $SORTED # TODO Is this how bamToFastq works?
		rm $SORTED

		printProgress "[extractFastqFromBAM] Compressing and moving extracted fastq file(s)..."
		for EX_FASTQ in ${OUTPUT//#.fastq/*.fastq}; do
			( gzip $EX_FASTQ
				mv $EX_FASTQ.gz ${EX_FASTQ//$TEMP/$CURRENT_DIRECTORY\/$FASTQ_DIRECTORY}.gz
			) &
			COMPRESS_PID_ARRAY[$COMPRESS_COUNTER]=$!
			((COMPRESS_COUNTER++))
		done
	done

	for pid in ${COMPRESS_PID_ARRAY[*]}
	do
		wait $pid
	done

	printProgress "[extractFastqFromBAM] Fastq extraction from $BAM_NAME completed."
	
	cd $CURRENT_DIRECTORY
}



########################################
#			 ALLELE-SPECIFIC FUNCTIONS			 # TODO move these into separate file
########################################

function checkPseudogenome () {
	if $ALLELE_SPECIFIC  ; then

		cd $TEMP_DIR
		
		CROSS_LIST=$(awk '(NF==4 && ($3!="" || $4!="")) {print $3"_"$4}' $INPUT_FILE | sort | uniq) #create array of crosses
		EXIT_SCRIPT=false

		for CROSS in $CROSS_LIST; do
			HAPLO_1=${CROSS//_*}
			HAPLO_2=${CROSS##*_}

			if [[ -z $HAPLO_1 || -z $HAPLO_2 ]]; then
				echo -e "ERROR:\t One of your samples only has 1 haplotype entered. \nPlease enter both haplotypes of allele pipeline" 
				EXIT_SCRIPT=true
				continue #continue to next iteration of cross
			fi

			MAKE_GENOME=true
			
			for DIP in $DIPLOID_GENOME_DIR/*; do #searching for folder with name of both haplotypes
				if [[ $DIP == *$HAPLO_1* && $DIP == *$HAPLO_2* ]]; then
					MAKE_GENOME=false

					if [[ ! -f $HAPLOID_GENOME_DIR/$HAPLO_1/$HAPLO_1".fa.refmap" ]] || \
						 [[ ! -f $HAPLOID_GENOME_DIR/$HAPLO_2/$HAPLO_2".fa.refmap" ]]; then #check RefMaps exist
						echo -e "ERROR:\t One of the refmap ($HAPLO_1 or $HAPLO_2) doesn't exists."
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
#	NAME=$SEARCH_KEY TODO this is probably broken now, will fix when allelic alignments needed
	printProgress "[setPseudogenome] Starting..."
	
	HAPLO_1=$(grep -e $NAME $INPUT_FILE | cut -f3 | head -n 1)
	HAPLO_2=$(grep -e $NAME $INPUT_FILE | cut -f4 | head -n 1)

	if [[ -z $HAPLO_1 || -z $HAPLO_2 ]]; then #if either haplotype is not entered, no allele specific for this data
		printProgress "[setPseudogenome] Haplotypes not provided, allele specific pipeline is not ran on $NAME"
		continue #dont continue the function for this SRACODE entry (move on to next code)
	fi
	
	CROSS_COMBO=$HAPLO_1"_"$HAPLO_2
	printProgress "[setPseudogenome] Haplotypes for $NAME are $HAPLO_1 and $HAPLO_2."
	
	for DIP in $DIPLOID_GENOME_DIR/*; do #searching for folder with name of both haplotypes
		if [[ $DIP == *$HAPLO_1* && $DIP == *$HAPLO_2* ]]; then
			GENOME_BUILD=$(basename $DIP)
			GENOME_FILE=$DIP/*".fa"
			STAR_GENOME_DIR=$DIP/*"-STAR"
			BOWTIE2_INDEXES=$DIP/$GENOME_BUILD

			HAPLO_1_REFMAP=$HAPLOID_GENOME_DIR/$HAPLO_1/$HAPLO_1".fa.refmap"
			HAPLO_2_REFMAP=$HAPLOID_GENOME_DIR/$HAPLO_2/$HAPLO_2".fa.refmap"

			printProgress "[setPseudogenome] Reference directory for $NAME is set to $DIP"
			break
		fi
	done
}

function unpackAllelic () { #working on bam that has already aligned to the pseudogenome (once per haplotype)
	cd $TEMP_DIR
	local HAPLO=$1
	
	printProgress "[unpackAllelic] Unpacking for $HAPLO"
	
	for TOT_RAW_BAM in $TEMP_DIR/*"_raw.bam"; do #should be in replicates (NAME1_rep1_raw.bam)
#take SAM file, align to pseudogenome
#split header into two, to separate the reads into two haplotypes, rename chr from hap1_chr to chr
		NAME_HAPLO_SUFFIX=$SEARCH_KEY"_"$HAPLO"_q$MIN_MAPQ"${TOT_RAW_BAM//$SEARCH_KEY/} #inserting haplotype name (-> NAME1_HAPLO_qMINMAPQ_rep1_raw.bam)
		NAME=${NAME_HAPLO//_raw.bam/} #will include _rep if applicable
		FILE_RAW_BAM=$NAME"_raw.bam"

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

		printProgress "[unpackAllelic $HAPLO] Finished unpacking $TOT_RAW_BAM for $HAPLO -> $FILE_RAW_BAM"
		
		refineBAM
#		cat $FILE"_markDupeMetrics.txt" >> $SEARCH_KEY"_alignLog.txt"

		rm $TOT_RAW_BAM
	done

	printProgress "[unpackAllelic] Unpacking bams for $HAPLO completed."
	cd $CURRENT_DIRECTORY
}

function projectAllelic () {
	if [[ $ALLELE_SPECIFIC == false ]]; then
		return 0 #dont run this code if the script is not allele specific
	fi 

	cd $TEMP_DIR

	printProgress "[projectAllelic] Starting..."
		
	for FILE_BAM in $CURRENT_DIRECTORY/$BAM_FOLDER_NAME/*".bam"; do
		if [[ $FILE_BAM != *$HAPLO_1* || $FILE_BAM != *$HAPLO_2* ]]; then #only projecting bams that were aligned to pseudogenome
			continue
		fi
		
		printProgress "[projectAllelic] Removing duplicates from $FILE_BAM with F=$FLAG"
		local NAME_MAP_FLAG=${FILE_BAM//.bam/}"_F"$FLAG
		FLAG_BAM=$NAME_MAP_FLAG".bam"
		$SAMTOOLS view -bh -F $FLAG $FILE_BAM > $FLAG_BAM
 
		#nonscaled projection
		printProgress "[projectAllelic] Converting $FLAG_BAM to bedGraph"
		$BEDTOOLS genomecov -ibam $FLAG_BAM -bg -split -scale $SCALING_FACTOR > $NAME_MAP_FLAG"_preProject.bedgraph"
		prepWigAndProject $NAME_MAP_FLAG $NAME_MAP_FLAG"_preProject.bedgraph"
			 
		#RPM scaled projection
		local READ_COUNT=$(samtools view -c $FLAG_NAME".bam")
		local SCALING_FACTOR=$(echo "scale=25; 1000000/$READ_COUNT" | bc) #calculating with 25 decimal places at least
		printProgress "[projectAllelic RPM] Detected $READ_COUNT filtered reads"

		if [[ $STRANDED_ALLELIC ]]; then #stranded RPM
			printProgress "[projectAllelic stranded RPM] Splitting $FLAG_BAM reads by 1st in pair..."
			local FIRST_PAIR_BAM=$NAME_MAP_FLAG"_firstInPair.bam"
			$SAMTOOLS view -bh -f 0x0040 $FLAG_BAM > $FIRST_PAIR_BAM		
			$BEDTOOLS genomecov -ibam $FIRST_PAIR_BAM -bg -split -scale $SCALING_FACTOR -strand + > $NAME_MAP_FLAG"_first_pos_RPM.bedGraph"
			$BEDTOOLS genomecov -ibam $FIRST_PAIR_BAM -bg -split -scale $SCALING_FACTOR -strand - > $NAME_MAP_FLAG"_first_neg_RPM.bedGraph"
			rm $FIRST_PAIR_BAM

			printProgress "[projectAllelic stranded RPM] Splitting $FLAG_BAM reads by 2nd in pair..."
			local SECOND_PAIR_BAM=$NAME_MAP_FLAG"_secondInPair.bam"
			$SAMTOOLS view -bh -f 0x0080 $FLAG_BAM > $SECOND_PAIR_BAM
			$BEDTOOLS genomecov -ibam $SECOND_PAIR_BAM -bg -split -scale $SCALING_FACTOR -strand + > $NAME_MAP_FLAG"_second_pos_RPM.bedGraph"
			$BEDTOOLS genomecov -ibam $SECOND_PAIR_BAM -bg -split -scale $SCALING_FACTOR -strand - > $NAME_MAP_FLAG"_second_neg_RPM.bedGraph"
			rm $SECOND_PAIR_BAM

			#plus/pos strand
			printProgress "[projectAllelic stranded RPM] Combining stranded bedgraphs for plus strand..."
			$BEDTOOLS unionbedg -i $NAME_MAP_FLAG"_first_neg_RPM.bedGraph" $NAME_MAP_FLAG"_second_pos_RPM.bedGraph" > $NAME_MAP_FLAG"_p_tmp.bedGraph"
			rm $NAME_MAP_FLAG"_first_neg_RPM.bedGraph" $NAME_MAP_FLAG"_second_pos_RPM.bedGraph"
			awk '{OFS="\t";FS="\t"} {print $1, $2, $3, $4+$5}' $NAME_MAP_FLAG"_p_tmp.bedGraph" > $NAME_MAP_FLAG"_pos_preProject.bedGraph"
			rm $NAME_MAP_FLAG"_p_tmp.bedGraph"
			prepWigAndProject $NAME_MAP_FLAG"_RPM_pos" $NAME_MAP_FLAG"_pos_preProject.bedGraph" " plus stranded RPM"

			#minus/neg strand
			printProgress "[projectAllelic stranded RPM] Combining stranded bedgraphs for minus strand..."
			$BEDTOOLS unionbedg -i $NAME_MAP_FLAG"_first_pos_RPM.bedGraph" $NAME_MAP_FLAG"_second_neg_RPM.bedGraph" > $NAME_MAP_FLAG"_n_tmp.bedGraph"
			rm $NAME_MAP_FLAG"_first_pos_RPM.bedGraph" $NAME_MAP_FLAG"_second_neg_RPM.bedGraph"
			awk '{OFS="\t";FS="\t"} {print $1, $2, $3, $4+$5}' $NAME_MAP_FLAG"_n_tmp.bedGraph" > $NAME_MAP_FLAG"_neg_preProject.bedGraph"
			rm $NAME_MAP_FLAG"_n_tmp.bedGraph"
			prepWigAndProject $NAME_MAP_FLAG"_RPM_neg" $NAME_MAP_FLAG"_pos_preProject.bedGraph" " minus stranded RPM"

		else #unstranded RPM
			printProgress "[projectAllelic unstranded RPM] Converting $FLAG_BAM to bedGraph"
			$BEDTOOLS genomecov -ibam $FLAG_BAM -bg -split -scale $SCALING_FACTOR > $NAME_MAP_FLAG"_RPM_preProject.bedgraph"
			prepWigAndProject $NAME_MAP_FLAG"_RPM" $NAME_MAP_FLAG"_RPM_preProject.bedgraph" " unstranded RPM"
		fi
	done

	printProgress "[projectAllelic] All haplotype-specific bams have completed projection."
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



########################################
#					TRACK-HUB FUNCTIONS				   #
########################################

function masterTrackHub () {



	mkdir -p $TRACK_HUB_DIR

	printf "hub <HubNameWithoutSpace>\nshortLabel <max 17 char, display on side>\nlongLabel Hub to display <fill> data at UCSC\ngenomesFile genomes.txt\nemail <email-optional>" > ./$TRACK_HUB_DIR/hub.txt

	TRACK_FOLDER=$TRACK_HUB_DIR/$GENOME_BUILD
	mkdir -p $TRACK_FOLDER
	TRACKDB=$TRACK_FOLDER/"trackDb.txt"
	printf "genome "$GENOME_BUILD"\ntrackDb "$GENOME_BUILD"/trackDb.txt" > $TRACK_HUB_DIR/genomes.txt		

	
	BAM_COVERAGE_ARGUMENTS="--binSize $BIN_SIZE -p $RUN_THREAD --normalizeUsing $NORMALIZE --smoothLength $SMOOTH_WIN --outFileFormat bigwig --minMappingQuality $MIN_MAPQ --ignoreDuplicates"

	PRINTED_DIR=""

	printProgress "[masterTrackHub] Starting..."
	
	for BAM_FILE in ./*/*.bam; do #currently in $CURRENT_DIRECTORY
		FOLDER_FILE=${BAM_FILE//.\//} #getting rid of the "./"
		FILE=$(basename $BAM_FILE) #leaving just the basename of the file
		FILE_NAME=${FILE//.bam/_}$NORMALIZE ##basename without file extension

		printProgress "[masterTrackHub] Bin size is $BIN_SIZE"
		if [[ $BIN_SIZE != 1 ]] ; then
			FILE_NAME=$FILE_NAME"_b"$BIN_SIZE
		fi

		printProgress "[masterTrackHub] Smoothing Window is $SMOOTH_WIN"
		if [[ $SMOOTH_WIN != 0 ]] ; then
			FILE_NAME=$FILE_NAME"_s"$SMOOTH_WIN
		fi

		FOLDER_NAME=${FOLDER_FILE%%\/*} #removing longest text of the matching pattern (after "/" in this case)

		printProgress "[masterTrackHub] Name of supertrack: $FOLDER_NAME"
		if [[ $FOLDER_NAME != $PRINTED_DIR ]] ; then # Create supertrack for housing associated tracks
			printf "track $FOLDER_NAME \nsuperTrack on show\nshortLabel $FOLDER_NAME \nlongLabel $FOLDER_NAME \nwindowingFunction mean\n\n" | tee -a $TRACKDB 
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
		elif [[ $FILE == *"BSSeq"* ]] || [[ $FILE == *"RRBS"* ]] || [[ $FILE == *"PBAT"* ]]; then
			FILE_NAME=${FILE//_$NORMALIZE/} ##remove unused normalization from name
			generateBSTrack
		else #ChIPseq - not stranded
			case $FILE_NAME in
				*K4me1*)
					COLOUR="0,100,255"
					;;
				*K4me3*)
					COLOUR="0,0,255"
					;;
				*K9me2*)
					COLOUR="200,100,100"
					;;
				*K9me3*)
					COLOUR="200,0,0"
					;;
				*K27me3*)
					COLOUR="255,0,150"
					;;
				*K27ac*)
					COLOUR="0,100,100"
					;;
				*K36me2*)
					COLOUR="100,0,100"
					;;
				*K36me3*)
					COLOUR="200,0,200"
					;;
				*H2AK119ub*)
					COLOUR="0,200,0"
					;;
				*PolII*)
					COLOUR="50,50,200"
					;;
				*)
					COLOUR="100,100,100"
					;;
			esac
			generateBigwigsUnstranded $FOLDER_FILE $FILE_NAME
			printTrackHubUnstranded $FOLDER_NAME $FILE_NAME $COLOUR
		fi
		
	done
}

function generateRNATrack () {
	COLOUR="50,0,200"
	if [[ $PAIRED == true ]] ; then
		echo "Extracting F reads over Actb..."
		$SAMTOOLS view -L $ACTB_BED -f 64 $FOLDER_FILE > Actb.sam
	else
		echo "Extracting reads over Actb..."
		$SAMTOOLS view -L $ACTB_BED $FOLDER_FILE > Actb.sam
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
		printTrackHubUnstranded $FOLDER_NAME $FILE_NAME $COLOUR
	else
		FILE_BIGWIG_POS=$TRACK_FOLDER/$FILE_NAME"_pos.bw"
		FILE_BIGWIG_NEG=$TRACK_FOLDER/$FILE_NAME"_neg.bw"
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
		printTrackHubStranded $FOLDER_NAME $FILE_NAME $COLOUR
	fi
}

function generateBSTrack () {
	COLOUR="0,0,0"
	FILE_TEMP_1=$TEMP_DIR/$FILE_NAME"_temp.bam"
	TEMP_BEDGRAPH=${FILE_TEMP_1//.bam/.bedGraph}
	TEMP_BEDGRAPH_2=${TEMP_BEDGRAPH//.bedGraph/_2.bedGraph}

	printProgress "[masterTrackHub generateBSTrack] Filtering $FILE for mapping quality of $MIN_MAPQ"
	$SAMTOOLS view -bh -@ $RUN_THREAD -q $MIN_MAPQ -o $FILE_TEMP_1 $FOLDER_FILE
	$BISMARK_METH_EXTRACT --gzip --multicore $BISMARK_THREAD --bedGraph --genome_folder $BISMARK_GENOME_DIR -o $TEMP_DIR $FILE_TEMP_1

	gunzip $TEMP_BEDGRAPH.gz
#		grep ^[^\*] $FILE_BEDGRAPH > $TEMP_DIR/temp.bedgraph
	tail -n +2 $TEMP_BEDGRAPH > $TEMP_BEDGRAPH_2
	sort -k1,1 -k2,2n $TEMP_BEDGRAPH_2 > $TEMP_BEDGRAPH
	rm $TEMP_BEDGRAPH_2
	$BEDGRAPHTOBW $TEMP_BEDGRAPH $CHROM_SIZES $TRACK_FOLDER/$FILE_NAME.bw
	rm $TEMP_BEDGRAPH
	printTrackHubUnstranded $FOLDER_NAME $FILE_NAME $COLOUR
	rm $FILE_TEMP_1
	rm $TEMP_DIR/CHG*
	rm $TEMP_DIR/CHH*
#	rm $TEMP_DIR/CpG* # Should we move this file?
	rm $TEMP_DIR/$FILE_NAME.*
#		rm $TEMP_DIR/*.bai
	#TODO Optional: keep more information? Lots of stuff being discarded
}

#$1=Name of filtered .bam file $2=Final name
function generateBigwigsUnstranded () {
	echo "Generating bigwig files..."
	$SAMTOOLS index $1
	$BAMCOVERAGE $BAM_COVERAGE_ARGUMENTS -b $1 --outFileName $TRACK_FOLDER/"$2.bw"
}

#$1=Supertrack name $2=Track name $3=Colour(RGB)
function printTrackHubUnstranded () { 
	printf "\ttrack %s\n\tparent %s\n\tshortLabel %s\n\tlongLabel %s\n\ttype bigWig\n\tbigDataUrl %s\n\tcolor %s\n\tvisibility full\n\tmaxHeightPixels 100:60:25\n\tautoScale on\n\talwaysZero on\n\n" $2 $1 $2 $2 $2".bw" $3 | tee -a $TRACKDB
}

#$1=Supertrack name $2=Track name $3=Colour(RGB)
function printTrackHubStranded () {
	printf "\ttrack %s\n\tparent %s\n\tcontainer multiWig\n\tshortLabel %s\n\tlongLabel %s\n\ttype bigWig\n\tvisibility full\n\tmaxHeightPixels 100:60:25\n\tconfigurable on\n\tautoScale on\n\talwaysZero on\n\taggregate transparentOverlay\n\tshowSubtrackColorOnUi on\n\tpriority 1.0\n\n" $2 $1 $2 $2 | tee -a $TRACKDB
	printf "\t\ttrack %s\n\t\tparent %s\n\t\tshortLabel %s\n\t\tlongLabel %s\n\t\ttype bigWig\n\t\tbigDataUrl %s\n\t\tcolor %s\n\t\tautoScale on\n\n" $2"_pos" $2 $2"_pos" $2"_pos" $2"_pos.bw" $3 | tee -a $TRACKDB
	printf "\t\ttrack %s\n\t\tparent %s\n\t\tshortLabel %s\n\t\tlongLabel %s\n\t\ttype bigWig\n\t\tbigDataUrl %s\n\t\tcolor %s\n\t\tautoScale on\n\n" $2"_neg" $2 $2"_neg" $2"_neg" $2"_neg.bw" $3 | tee -a $TRACKDB
}

