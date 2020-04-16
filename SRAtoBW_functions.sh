#! /bin/bash

###List of functions that can be used in various pipelines for download, align and creating track hub

##System Specific Configuration

#TODO: Only if you ever plan on using these functions directly in the command line (otherwise ignore)
#Ensure correct config file if using
if [[ -z $SCRIPTS_DIR ]]; then
  pushd $(dirname $0) > /dev/null
  SCRIPTS_DIR=$(pwd -P)
  popd > /dev/null
  source $SCRIPTS_DIR/BRC.config
fi


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
SEP_PARA=false
STRANDED_ALLELIC=false
TRIM_READ=false
USE_BOWTIE=false
USE_SERVER=false

CODE_ARRAY=""
BIN_SIZE=1
FLAG=1540 #read unmapped, read fails platform/vendor quality checks, read is PCR or optical duplicate 
GENOME_BUILD="mm10"
MIN_MAPQ=5
NORMALIZE="CPM"
SMOOTH_WIN=0

DEPENDENCIES=($ESEARCH $EFETCH $FASTERQDUMP "$TRIMMOMATIC" $STAR $BISMARK $BOWTIE2 $BWA $SAMTOOLS "$PICARD" awk $BAM2FASTQ $BEDGRAPHTOBW $BAMCOVERAGE)

# Help Menu
OPTIONS="hi:ab:B:d:Df:Fg:km:M:n:N:ors:St:Tux:X"

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
-m\tMemory to give each thread (Format=XG/M).\n\t
-M\tMinimum mapping quality for bigwig generation. Default=5\n\t
-n\tBin size for bigwig generation. Larger bins to smooth noisy\n\t\tdata. Default=1\n\t
-N\tNormalization method for bigwigs. Accepted: CPM, RPKM\n\t\t(Default=CPM)\n\t
-r\tKeep replicates after collapsing. Default=false.\n\t
-s\tObtained stranded RPM tracks for allele-specific runs.\n\t\tDefault=false\n\t
-S\tUse of server to submit parallel jobs. Default=FALSE\n\t
-t\tNumber of Threads to use. Default=6 (Check in config file)\n\t
-T\tTrim .fastq files after download.\n\t
-u\tChanging ChIPseq aligner to bowtie2. Default=BWA\n\t
-w\tSmoothing window. Will smooth bigwigs in a rolling window of\n\t\tthis size. Default=0\n\t"



########################################
#               FUNCTIONS              #
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
      f) #use exisiting fastq files for downstream modules in pipeline (will skip download)
        FASTQ_INPUT=true
        if $FASTQ_ONLY ; then
          echo -e "ERROR:\tFiles already in Fastq format."
          exit 1
        fi
        FASTQ_DIRECTORY=${OPTARG}
        ;;
      F) #Stop pipeline after downloading
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
      S)
        USE_SERVER=true
        ;;
      t)
        RUN_THREAD=${OPTARG}
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
      x)
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

function setUp () { #set up log file for parallel runs
  if [[ $SEP_PARA == false ]]; then 
    LOG_FILE=$CURRENT_DIRECTORY/$SEARCH_KEY"_"$(date '+%y-%m-%d')"_log.txt"
    checkDependencies
    printProgress "[setUp] Script started at [$(date)]"
    printProgress "[setUp] Search key for the set: $SEARCH_KEY"
    printProgress "[setUp] SRA array: $CODE_ARRAY"
    printProgress "[setUp setGenome] Genome used for data is $GENOME_BUILD"
  fi
}

function printProgress () {
  echo -e $1 | tee -a $LOG_FILE #using tee will also show in stdout
}

function checkFileExists () {
  if [[ ! -f $1 ]]; then
    echo "ERROR:\tFile $1 does not exist"
    exit 1
  fi
}

### Create neccessary files for reference genome
function setGenome () {
	if [[ $FASTQ_ONLY == false ]] ; then
		mkdir -p $TRACK_HUB_DIR
		printf "hub <HubNameWithoutSpace>\nshortLabel <max 17 char, display on side>\nlongLabel Hub to display <fill> data at UCSC\ngenomesFile genomes.txt\nemail <email-optional>" > ./$TRACK_HUB_DIR/hub.txt
	fi
  
  MOUSE="." #TODO
  RAT="."   #CURRENTLY BRC LOCATION!

  case $1 in
		"mm9")
			GENOME_DIR=$GENOME_DIR/$MOUSE/$1/
			BISMARK_GENOME_DIR=$GENOME_DIR
			BOWTIE_GENOME_DIR=$GENOME_DIR
			if [[ $FASTQ_ONLY == false ]]; then
				printf "chr5\t143666090\t143666091" > Actb.bed
				TRACK_FOLDER=$TRACK_HUB_DIR/$1
				mkdir -p $TRACK_FOLDER
				TRACKDB=$TRACK_FOLDER/"trackDb.txt"
				printf "genome mm9\ntrackDb mm9/trackDb.txt" > $TRACK_HUB_DIR/genomes.txt
			fi
			;;
		"mm10")
		  GENOME_DIR=$GENOME_DIR/$MOUSE/$1/ 
			if [[ $FASTQ_ONLY == false ]] ; then
				printf "chr5\t142904365\t142904366" > Actb.bed
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
		"rn5")
		  GENOME_DIR=$GENOME_DIR/$RAT/$1/
			if [[ $FASTQ_ONLY == false ]] ; then
				printf "chr12\t15748011\t15748012" > Actb.bed
				TRACK_FOLDER=$TRACK_HUB_DIR/$1
				mkdir -p $TRACK_FOLDER
				TRACKDB=$TRACK_FOLDER/"trackDb.txt"
				printf "genome rn5\ntrackDb rn5/trackDb.txt" > $TRACK_HUB_DIR/genomes.txt
			fi
			;;
    "hg19")
      GENOME_DIR=$GENOME_DIR/"Hsa"/$1/
			if [[ $FASTQ_ONLY == false ]] ; then
				printf "chr7\t5527531\t5527532" > Actb.bed
				TRACK_FOLDER=$TRACK_HUB_DIR/$1
				mkdir -p $TRACK_FOLDER
				TRACKDB=$TRACK_FOLDER/"trackDb.txt"
				printf "genome hg19\ntrackDb hg19/trackDb.txt" > $TRACK_HUB_DIR/genomes.txt
			fi
			;;
		"oryCun2")
		  GENOME_DIR=$GENOME_DIR/"oryCun"/$1/
			if [[ $FASTQ_ONLY == false ]] ; then
				printf "chr7\t87232876\t87232877" > Actb.bed
				TRACK_FOLDER=$TRACK_HUB_DIR/$1
				mkdir -p $TRACK_FOLDER
				TRACKDB=$TRACK_FOLDER/"trackDb.txt"
				printf "genome oryCun2\ntrackDb oryCun2/trackDb.txt" > $TRACK_HUB_DIR/genomes.txt
			fi
			;;
		"mesAur1")
		  GENOME_DIR=$GENOME_DIR/"mesAur"/$1/
			if [[ $FASTQ_ONLY == false ]] ; then
				printf "KB708222.1\t2764075\t2764076" > Actb.bed
				TRACK_FOLDER=$TRACK_HUB_DIR/$1
				mkdir -p $TRACK_FOLDER
				TRACKDB=$TRACK_FOLDER/"trackDb.txt"
				printf "genome mesAur1\ntrackDb mesAur1/trackDb.txt" > $TRACK_HUB_DIR/genomes.txt
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

### Create subset of SRACODE array to be used in parallel running
function parallelRun () {
  if $SEP_PARA; then
    echo -e "Separating input file into subsets for parallel runs..."

    INPUT_FILE=$CURRENT_DIRECTORY/$FILENAME
    createCodeArray
    CURRENT_SET=""
  
    for code in $CODE_ARRAY; do
      SRACODE=$code
      NAME=$(grep -e $SRACODE $INPUT_FILE | cut -f2)

      if [[ $NAME == *_[Rr]ep* ]] ; then        
        if [[ $CURRENT_SET != ${NAME//_[Rr]ep*/}* ]]; then
          CURRENT_SET=${NAME//_[Rr]ep*/}
          declare -a SUB_ARRAY=$(grep -e $CURRENT_SET $INPUT_FILE | cut -f1)
          echo "calling "$(basename $SHELL_SCRIPT) "on" $CURRENT_SET
          
          if $USE_SERVER; then
            echo "Submitting on server"
            DATE=$(date '+%y-%m-%d')
            $SERVER_SUBMIT "MasterDAT_"$DATE"_"$CURRENT_SET $SHELL_SCRIPT $PASS_ARG -x $CURRENT_SET -X $SUB_ARRAY
            sleep 30 #pause for 30 secs before running next code b/c fetching data takes some time
          else
            $SHELL_SCRIPT $PASS_ARG -x $CURRENT_SET -X $SUB_ARRAY & 
            wait $! #wait for the script above to finish running before moving onto the next set (avoid overload)
          fi

        fi

      else
        declare -a SUB_ARRAY=$SRACODE
        echo "calling "$(basename $SHELL_SCRIPT) "on" $NAME

        if $USE_SERVER ; then
          echo "Submitting on server"
          DATE=$(date '+%y-%m-%d')
          $SERVER_SUBMIT "MasterDAT_"$DATE"_"$NAME $SHELL_SCRIPT $PASS_ARG -x $NAME -X $SUB_ARRAY
          sleep 30
        else
          $SHELL_SCRIPT $PASS_ARG -x $NAME -X $SUB_ARRAY
          wait $!
        fi

      fi
    done
    exit #exit the script b/c don't want to run the rest of the code on every single SRACODE again
  fi

}

### Downloading files specified from the tab-delimited file to fastq files
### When calling function by itself $1=SRA input file
function masterDownload () {
  if $BAM_INPUT || $FASTQ_INPUT; then #if input is either BAM/FASTQ then no download is needed - exit masterDownload function
    return 0
  fi

  mkdir $FASTQ_DIRECTORY
  if $BAM_REALIGNMENT; then
    obtainFastqFromBAM
    return 0 #after extracting fastq from BAM, can exit masterDownload function
  fi

  INPUT_FILE=$CURRENT_DIRECTORY/$FILENAME
  
  if [[ -z $CODE_ARRAY ]]; then #this will basically only be used if the function was called by itself (may delete later)
    INPUT_FILE=$CURRENT_DIRECTORY/$F1
    createCodeArray
  fi
  
  cd $TEMP_DIR  

  printProgress "[masterDownload] Started at [$(date)]"
    
  for code in $CODE_ARRAY; do
    downloadReads
    extractType
    extractPaired
    extractFastq
  done
  cd $CURRENT_DIRECTORY

  printProgress "[masterDownload] All fastq files for $SEARCH_KEY is downloaded at [$(date)]"

  if $FASTQ_ONLY; then
    printProgress "[masterDownload] Fastq files only. Exit script."
    exit 0 #exit the whole script b/c only requires fastq files
  fi
}
 
### Create an array of SRACODE passed listed in the tab-delimited file
function createCodeArray () {
  while read n; do #read command will read a line of and split them into words (2 words in files)
  	declare -a m=($n) #-a for an array m to be stored based on the line (SRACODE, NAME)
  	CODE_ARRAY=$CODE_ARRAY" "${m[0]} #add only SRACODE to this list
  done < $INPUT_FILE #feeding FILE as stdin to this while read
}


### Download reads from SRA
function downloadReads () {
  DL_COUNTER=0

  SRACODE=$code
  NAME=$(grep -e $SRACODE $INPUT_FILE | cut -f2)

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
  
  printProgress "[masterDownload] Finished $SRACODE -> $NAME reads download."
}

### Determine data type (could be called with just SRACODE)
function extractType() {
  SRACODE=${1:-$SRACODE} #giving option to just use function to see what type of sequence it is? is this neccessary?

  TYPE=$($ESEARCH -db sra -query $SRACODE \
        | $EFETCH --format runinfo \
        | cut -d ',' -f 13 \
        | tail -n 2 \
        | head -n 1)
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

### Extract fastq using fasterq-dump, compressed then concatenate fastq
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
	  	cat $SEARCH_DL*_1.fastq.gz > $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_1.fastq.gz"
	  	cat $SEARCH_DL*_2.fastq.gz > $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_2.fastq.gz"
  	else # If only single file, move, don't copy
	  	mv $SEARCH_DL*_1.fastq.gz $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_1.fastq.gz"
	  	mv $SEARCH_DL*_2.fastq.gz $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_2.fastq.gz"
	  fi
 
	else #Single-End Reads
    COUNT=$(ls -1 $SEARCH_DL*.fastq.gz | wc -l)
    if [[ $COUNT > 1 ]] ; then # If only single file, move, don't copy
      cat $SEARCH_DL*.fastq.gz > $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME".fastq.gz"
    else
      mv $SEARCH_DL*.fastq.gz $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME".fastq.gz"
    fi
	fi

	printProgress "[masterDownload] Fastq files for $NAME are moved to $CURRENT_DIRECTORY/$FASTQ_DRECTORY"
  rm $SEARCH_DL
}

### Trim fastq files for adaptors using trimmomatic
### Could be called with $1=directory that holds fastq files
function trimReads () {
  if [[ $TRIM_READ == false || -z $1 ]]; then #exit script if not needed
    printProgress "[trimReads] Reads not trimmed."
    return 0
  fi

  cd $TEMP_DIR
  printProgress "[trimReads] Started at [$(date)]"
    
  for FILE in $CURRENT_DIRECTORY/${1:-$FASTQ_DIRECTORY}/*$SEARCH_KEY*fastq.gz; do
    determinePairedFastq

    if [[ $PAIRED_END ]]; then
      printProgress "[trimReads] Trimming "$NAME"_1.fastq.gz and "$NAME"_2.fastq.gz..."

      $TRIMMOMATIC PE -threads 10 $FASTQ_PATH"_1.fastq.gz" $FASTQ_PATH"_2.fastq.gz" \
      $NAME"_trim_1.fastq.gz" $NAME"_unpaired_trim_1.fastq.gz" $NAME"_trim_2.fastq.gz" $NAME"_unpaired_trim_2.fastq.gz" \
      ILLUMINACLIP:$ILLUMINA_ADAPATORS_ALL":2:30:10" \
      SLIDINGWINDOW:4:20 \
      MINLEN:36        
        
      mv $NAME"_trim_1.fastq.gz" $FASTQ_PATH"_1.fastq.gz"
      mv $NAME"_trim_2.fastq.gz" $FASTQ_PATH"_2.fastq.gz"

    else #Single-End
      printProgress "[trimReads] Trimming $NAME.fastq.gz"

      $TRIMMOMATIC SE -threads 10 $FASTQ_PATH".fastq.gz" $FASTQ_PATH"_trim.fastq.gz" \
      ILLUMINACLIP:$ILLUMINA_ADAPATORS_ALL":2:30:10" \
      SLIDINGWINDOW:4:20 \
      MINLEN:36

		  mv $NAME"_trim.fastq.gz" $FASTQ_PATH".fastq.gz"
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

  else
    FASTQ_PATH=${FILE//.fastq.gz}
    NAME=${FASTQ_PATH##*/}

    PAIRED_END=false
    FILE_FASTQ=$FASTQ_PATH".fastq.gz"
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
### $1 can be a specific dataset prefix, will only align according to that tag
function masterAlign () {
  if $BAM_INPUT; then #if input is BAM alignment is not needed - exit function
    return 0
  fi

  cd $TEMP_DIR
  printProgress "[masterAlign] Started at [$(date)]"
  
  STAR_ARGUMENTS="--genomeDir $STAR_GENOME_DIR --runThreadN $RUN_THREAD --sjdbOverhang 70 --outFilterType BySJout --twopassMode Basic --twopass1readsN 1000000000 --outSAMunmapped Within --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --readFilesCommand zcat "

  BISMARK_ARGUMENTS="--chunkmbs $BISMARK_MEM --multicore $BOWTIE_THREAD --genome $BISMARK_GENOME_DIR "
  SEARCH_KEY=${1:-$SEARCH_KEY}

  for FILE in $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/*$SEARCH_KEY*fastq.gz; do
    determinePairedFastq

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

    if $ALLELE_RUN; then #allele specific run will need to do unpacking
      continue #move onto next set instead of refining the produced RAW BAM
    fi
    
    refineBam
    
  done

  printProgress "[masterAlign] Alignment of $SEARCH_KEY fastq files to $GENOME_BUILD completed at [$(date)]"

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

  printProgress "[masterAlign STAR] Alignment completed -> $FILE_RAW_BAM"
  mv $FILE_STAR_OUTPUT $FILE_RAW_BAM
  rm -r $SEARCH_KEY*"STAR"*
}

### Alignment of BSSeq data using Bismark
function alignBismark () {
  BISMARK_OUTPUT=$NAME"_bismark_bt2.bam"

  if [[ $NAME == *"PBAT"* ]] ; then
    BISMARK_ARGUMENTS=$BISMARK_ARGUMENTS"--pbat "
  else
    BISMARK_ARGUMENTS=$BISMARK_ARGUMENTS"--non_directional "
  fi

  if $PAIRED_END; then
    printProgress "[masterAlign Bismark] Aligning "$NAME"_1.fastq.gz and "$NAME"_2.fastq.gz to genome..."
    $BISMARK $BISMARK_ARGUMENTS -1 $FILE_FASTQ1 -2 $FILE_FASTQ2
    BISMARK_OUTPUT=${BISMARK_OUTPUT//_bismark_bt2.bam/_1_bismark_bt2_pe.bam}

  else #Single-End
    printProgress "[masterAlign Bismark] Aligning $NAME.fastq.gz to genome..."
    $BISMARK $BISMARK_ARGUMENTS $FILE_FASTQ
  fi

  printProgress "[masterAlign Bismark] Alignment completed -> $FILE_RAW_BAM"
  mv $BISMARK_OUTPUT $FILE_RAW_BAM
  cleanBismark
  rm *report.txt
}

### Alignment of HiC data with HiCUP, assumed to be paired end 
function alignHiCUP () {
  ENZYME=$(grep -e $NAME $INPUT_FILE | cut -f5 | head -n 1) 

  if [[ -z $ENZYME ]]; then
    printProgress "[masterAlign HiCUP] ERROR:\tDigestion enzyme was not entered for $SEARCH_KEY HiC dataset..."
    exit 1
  fi
    
  printProgress "[masterAlign HiCUP] Digesting $GENOME_BUILD with $ENZYME"
  hicup_digester --re1 $ENYZME --outdir $TEMP_DIR  --genome $GENOME_BUILD"_"$SEARCH_KEY $GENOME_FILE
  local DIGEST_FILE=$(ls $TEMPDIR/Digest*$SEARCH_KEY*)

  printProgress "[masterAlign HiCUP] Aligning $NAME to $GENOME_BUILD"
  $HICUP --bowtie2 $BOWTIE2 --digest $DIGEST_FILE --index $GENOME_DIR/$GENOME_BUILD --outdir $BAM_FOLDER --temp $TEMP_DIR --threads $RUN_THREAD $NAME"_1.fastq.gz" $NAME"_2.fastq.gz"

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
  $SAMTOOLS sort -@ $RUN_THREAD -m $THREAD_MEM -o $FILE_SORTED_BAM -T $NAME $FILE_CLEANED_BAM
    
  printProgress "[refineBAM] Marking duplicates..." #not removing the duplicates
  $PICARD MarkDuplicates I=$FILE_SORTED_BAM O=$FILE_BAM M=$NAME"_markDupeMetrics.txt"

  mv $FILE_BAM $BAM_FOLDER
  printProgress "[refineBAM] Final $FILE_BAM is moved to $BAM_FOLDER."
  
  rm $FILE_RAW_BAM $FILE_CLEANED_BAM $FILE_SORTED_BAM #remove all the buffer bam files

}

### Searches for .bam files containing _Rep#.bam in nested directories and combines them into a single bam
function collapseReplicates () {
  cd $TEMP_DIR
	CURRENT=""
	printProgress "[collapseReplicates] Started at [$(date)]"
	
	for FILE in $CURRENT_DIRECTORY/*/*$SEARCH_KEY*.bam; do
		if [[ $FILE == *_[Rr]ep* ]]; then 
		  MERGED_BAM=${FILE//_[Rr]ep*.bam/.bam}
			if [[ $CURRENT != ${FILE//_[Rr]ep*.bam/}*.bam ]]; then 
				CURRENT=$FILE
				printProgress "[collapseReplicates] Merging to $MERGED_BAM"
				$SAMTOOLS merge -@ $RUN_THREAD $MERGED_BAM ${FILE//_[Rr]ep*.bam/}*.bam
				
        if [[ $KEEP_REPLICATES == false ]]; then
          printProgress "[collapseReplicates] Removing replicates"
				  rm ${FILE//_[Rr]ep*.bam/}_*ep*.bam
        fi

        printProgress "[collapseReplicates] Indexing BAM file..."
        $SAMTOOLS index ${MERGED_BAM//.bam/}*.bam #index merged & replicates (if available)

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

#    printProgress "[collapseReplicates] Obtaining flagstats for $MERGED_BAM flagstats"
#    $SAMTOOLS flagstat $MERGED_BAM
		
	done

  printProgress "[collapseReplicates] All replicates for $SEARCH_KEY has been combined at [$(date)]"

  cd $CURRENT_DIRECTORY
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
    NAME=${BAM_INPUT##*/}
    SORTED="temp_sorted.bam"
    mkdir -p "temp"
    OUTPUT="temp"/${NAME//.bam/#.fastq} #the # will be replaced by _1/_2 for PE reads in bam2fastq

    echo "Sorting $NAME by read name..."
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



########################################
#       ALLELE-SPECIFIC FUNCTIONS      #
########################################

function checkPseudogenome () {
  if $ALLELE_SPECIFIC && $SEP_PARA ; then

    cd $TEMP_DIR
    
    CROSS_LIST=$(awk '(NF==4 && ($3!="" || $4!="")) {print $3"_"$4}' $INPUT_FILE | sort | uniq) #create array of crosses
    EXIT_SCRIPT=false

    for CROSS in $CROSS_LIST; do
      HAPLO_1=${CROSS//_*}
      HAPLO_2=${CROSS##*_}

      if [[ -z $HAPLO_1 || -z $HAPLO_2 ]]; then
        echo -e "ERROR:\t One of your samples only have 1 haplotype entered. \nPlease ensure to enter both haplotypes of allele pipeline" \
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
  NAME=$SEARCH_KEY
  printProgress "[setPseudogenome] Started at [$(date)]"
  
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
  
  printProgress "[unpackAllelic] Unpacking for $HAPLO started at [$(date)]"
  
  for TOT_RAW_BAM in $TEMP_DIR/$SEARCH_KEY*"_raw.bam"; do #should be in replicates (NAME1_rep1_raw.bam)
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
#    cat $FILE"_markDupeMetrics.txt" >> $SEARCH_KEY"_alignLog.txt"

    rm $TOT_RAW_BAM
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



########################################
#          TRACK-HUB FUNCTIONS         #
########################################

function masterTrackHub () {
  BAM_COVERAGE_ARGUMENTS="--binSize $BIN_SIZE -p $RUN_THREAD --normalizeUsing $NORMALIZE --smoothLength $SMOOTH_WIN --outFileFormat bigwig --minMappingQuality $MIN_MAPQ --ignoreDuplicates"

  PRINTED_DIR=""
  
  for BAM_FILE in ./*/*$SEARCH_KEY*.bam; do #currently in $CURRENT_DIRECTORY
    FOLDER_FILE=${BAM_FILE//.\//} #getting rid of the "./"
    FILE=$(basename $BAM_FILE) #leaving just the basename of the file
    FILE_NAME=${FILE//.bam/_}$NORMALIZE ##basename without file extension

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
	  elif [[ $FILE == *"BSSeq"* ]] || [[ $FILE == *"RRBS"* ]] || [[ $FILE == *"PBAT"* ]]; then
      generateBSTrack
	  else #ChIPseq - not stranded
		  generateBigwigsUnstranded $FOLDER_FILE $FILE_NAME
		  printTrackHubUnstranded $FOLDER_NAME $FILE_NAME
	  fi
	  
  done
  rm Actb.bed
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

function generateBSTrack () {
# FILE_BEDGRAPH=$TEMP_DIR/$FILE_NAME".bedGraph"
	FILE_TEMP_1=$TEMP_DIR"/temp.bam"
	TEMP_BEDGRAPH=${FILE_TEMP_1//.bam/.bedGraph}
#		FILE_TEMP_2=$TEMP_DIR"/temp2.bam"
  echo "Filtering" $FILE "for mapping quality..."
  $SAMTOOLS view -bh -@ $RUN_THREAD -q $MIN_MAPQ -o $FILE_TEMP_1 $FOLDER_FILE
  $BISMARK_METH_EXTRACT --gzip --multicore $RUN_THREAD --bedGraph --genome_folder $BISMARK_GENOME_DIR -o $TEMP_DIR $FILE_TEMP_1
#		mv temp2.bedGraph.gz $FILE_BEDGRAPH.gz # TODO probably not correct
  gunzip $TEMP_BEDGRAPH.gz
#		grep ^[^\*] $FILE_BEDGRAPH > $TEMP_DIR/temp.bedgraph
#		mv $TEMP_BEDGRAPH $FILE_BEDGRAPH
#		sort -k1,1 -k2,2n $TEMP_DIR/temp.bedgraph > $FILE_BEDGRAPH
#		rm $TEMP_DIR/temp.bedgraph
  $BEDGRAPHTOBW $TEMP_BEDGRAPH $CHROM_SIZES $TRACK_FOLDER/$FILE_NAME.bw
  rm $TEMP_BEDGRAPH
  printTrackHubUnstranded $FOLDER_NAME $FILE_NAME
  rm $FILE_TEMP_1
  rm $TEMP_DIR/CHG*
  rm $TEMP_DIR/CHH*
  rm $TEMP_DIR/CpG*
  rm $TEMP_DIR/$FILE_NAME.*
#		rm $TEMP_DIR/*.bai
  #TODO Optional: keep more information? Lots of stuff being discarded
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
