#! /bin/bash

###List of functions that can be used in various pipelines for download, align and creating track hub

##System Specific Configuration

#TODO: alter any system specific variables and tools path through config file
#Ensure function and config files are within same directory (sourcing will not work otherwise)
source $FUNCTIONS_DIR/BRC.config

shopt -s nocaseglob #turning off case match for files 

## Non System Specific Variables
CURRENT_DIRECTORY=$(pwd)
FASTQ_DIRECTORY="./Fastq"

ALLELE_SPECIFIC=false
BAM_INPUT=false
BAM_REALIGNMENT=false
FASTQ_INPUT=false
FASTQ_ONLY=false
KEEP_FASTQ=false
KEEP_REPLICATES=false
PARALLEL=false
TRIM_READ=false
USE_BOWTIE=false
USE_SERVER=false

CODE_ARRAY=""
GENOME_BUILD="mm10"

BIN_SIZE=1
SMOOTH=0
NORMALIZE="CPM"
MIN_MAPPING_QUAL=5

DEPENDENCIES=($ESEARCH $EFETCH $FASTERQDUMP "$JAVA $TRIMMOMATIC" $STAR $BISMARK $BWA $SAMTOOLS "$JAVA $MARKDUPS" awk $BAM2FASTQ $BEDGRAPHTOBW $BAMCOVERAGE)

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
-m\tMemory to give each thread (Format=XG/M). Default=1G\n\t
-M\tMinimum mapping quality for bigwig generation. Default=5\n\t
-n\tBin size for bigwig generation. Larger bins to smooth noisy\n\t\tdata. Default=1\n\t
-N\tNormalization method for bigwigs. Accepted: CPM, RPKM\n\t\t(Default=CPM)\n\t
-r\tKeep replicates after collapsing. Default=false.\n\t
-s\tSmoothing window. Will smooth bigwigs in a rolling window of\n\t\tthis size. Default=0\n\t
-S\tUse of server to submit parallel jobs. Default=FALSE\n\t
-t\tNumber of Threads to use. Default=10\n\t
-T\tTrim .fastq files after download.\n\t
-u\tChanging ChIPseq aligner to bowtie2. Default=BWA\n\t"





############### FUNCTIONS ###############
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
        PARALLEL=true
        PASS_ARG=$@
        FILENAME=${OPTARG}
        ;;
      a)
        ALLELE_SPECIFIC=true
        ;;
      b) #bam input for trackhub
        BAM_INPUT=true
        if $FASTQ_ONLY ; then
          echo "Incompatible options. Cannot generate .fastq files from .bam files."
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
          echo "Files already in Fastq format"
          exit 1
        fi
        FASTQ_DIRECTORY=${OPTARG}
        ;;
      F) #Stop pipeline after downloading
        FASTQ_ONLY=true
        if $BAM_INPUT ; then
          echo "Incompatible options. Cannot generate .fastq files from .bam files."
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
        MIN_MAPPING_QUAL=${OPTARG}
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
        SMOOTH=${OPTARG}
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
      x)
        PARALLEL=false
        SEARCH_KEY=${OPTARG}
        echo "Search key for the set: $SEARCH_KEY"
        ;;
      X)
        CODE_ARRAY=$(echo $@ | sed 's/.*-X //')
        echo "SRA array: $CODE_ARRAY"
        ;;
      \?)
        echo -e "\n###############\nERROR: Invalid Option! \nTry '$(basename $0) -h' for help.\n###############" >&2
        exit 1
        ;;
    esac
  done
  setGenome $GENOME_BUILD  
  checkDependencies
  shopt -s nocasematch #turn off case matching
}

### Create subset of SRACODE array to be used in parallel running
function parallelRun () {
  if $PARALLEL; then

    INPUT_FILE=$CURRENT_DIRECTORY/$FILENAME
    createCodeArray
    CURRENT_SET=""
  
    for code in $CODE_ARRAY; do
      SRACODE=$code
      NAME=$(grep -e $SRACODE $INPUT_FILE | cut -f2)

      if [[ $NAME == *_Rep* ]] ; then
        CASE_REP=${NAME##*_} ##only want the part that say "[Rr]ep" for renaming
        CASE_REP="_"${CASE_REP//ep*/ep} #changing [Rr]ep# to _[Rr]ep
        
        if [[ $CURRENT_SET != ${NAME//_Rep*/$CASE_REP}* ]]; then
          CURRENT_SET=${NAME//_Rep*/$CASE_REP}*
          declare -a SUB_ARRAY=$(grep -e ${NAME//_Rep*/$CASE_REP}* $INPUT_FILE | cut -f1)
          echo "calling "$(basename $SHELL_SCRIPT) "on" ${CURRENT_SET//_Rep*/$CASE_REP}
          
          if $USE_SERVER; then
            echo "Submitting on server"
            $SERVER_SUBMIT "MasterDAT_"${CURRENT_SET//_Rep*/} $SHELL_SCRIPT $PASS_ARG -x ${CURRENT_SET//_Rep*/$CASE_REP} -X $SUB_ARRAY
            sleep 30 #pause for 30 secs before running next code b/c fetching data takes some time
          else
            $SHELL_SCRIPT $PASS_ARG -x ${CURRENT_SET//_Rep*/$CASE_REP} -X $SUB_ARRAY & 
            wait $! #wait for the script above to finish running before moving onto the next set (avoid overload)
          fi

        fi

      else
        declare -a SUB_ARRAY=$SRACODE
        echo "calling "$(basename $SHELL_SCRIPT) "on" $NAME

        if $USE_SERVER ; then
          echo "Submitting on server"
          $SERVER_SUBMIT "MasterDAT_"$NAME $SHELL_SCRIPT $PASS_ARG -x $NAME -X $SUB_ARRAY
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
 
  if [[ -z $CODE_ARRAY ]]; then #this will basically only be used if the function was called by itself (may delete later)
    INPUT_FILE=$CURRENT_DIRECTORY/$1
    createCodeArray
  fi
  
  cd $TEMP_DIR  

  for code in $CODE_ARRAY; do
    downloadReads
    extractType
    extractPaired
    extractFastq
  done
  cd $CURRENT_DIRECTORY

  if $FASTQ_ONLY; then
    echo "Fastq files only"
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

  echo "downloading $SRACODE to $TEMP_DIR..."
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
  
  echo "finished $SRACODE -> $NAME download!"
}

### Determine data type (could be called with just SRACODE)
function extractType() {
  SRACODE=${1:-$SRACODE} #giving option to just use function to see what type of sequence it is? is this neccessary?

  
  TYPE=$($ESEARCH -db sra -query $SRACODE \
        | $EFETCH --format runinfo \
        | cut -d ',' -f 13 \
        | tail -n 2 \
        | head -n 1)
  echo "Data type: "$TYPE
  
  if [[ $TYPE == "ChIP-Seq" ]] && \
     [[ $NAME != *"ChIP"* ]]; then
     NAME="ChIPseq_"$NAME

  elif [[ $TYPE == "RNA-Seq" ]] && \
       [[ $NAME != *"RNA"* ]] ; then
       NAME="RNAseq_"$NAME

  elif [[ $TYPE == "Bisulfite-Seq" ]] && \
       [[ $NAME != *"RRBS"* ]] && \
       [[ $NAME != *"BSSeq"* ]] && \
       [[ $NAME != *"PBAT"* ]] ; then
       NAME="BSSeq_"$NAME

  elif [[ $TYPE == "DNase-Hypersensitivity" ]] && \
       [[ $NAME != *"DNase"* ]]; then
       NAME="DNase_"$NAME
  fi
  echo "Name of data: "$NAME
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
    echo "Data are single-end."
  else
    PAIRED_END=true
    echo "Data are paired-end."
  fi
} 

### Extract fastq using fasterq-dump, compressed then concatenate fastq
function extractFastq () { 
  COMPRESS_COUNTER=0
  COMPRESS_PID_ARRAY=""

  local SEARCH_DL=*$SEARCH_KEY*[DSE]RR*

  echo "Dumping .fastq files..."
	for SRA_FILE in $SEARCH_DL; do
		$FASTERQDUMP -e $RUN_THREAD --split-files ./$SRA_FILE
	done

  echo "Compressing fastq files..." #simultaneously zip all the dump fastq files for that read
  for DL_FASTQ in $SEARCH_DL*fastq; do 
    gzip $DL_FASTQ & 
    COMPRESS_PID_ARRAY[$COMPRESS_COUNTER]=$!
    ((COMPRESS_COUNTER++))
  done
  
  for gzipid in ${COMPRESS_PID_ARRAY[*]}; do
    wait $gzipid #wait for all gzip to finish
  done

  if $PAIRED_END; then
  	echo "Concatenating .fastq files..."
	  COUNT=$(ls -1 $SEARCH_DL*_1.fastq.gz | wc -l)
  	if [[ $COUNT > 1 ]] ; then 
	  	cat $SEARCH_DL*_1.fastq.gz > $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_1.fastq.gz"
	  	cat $SEARCH_DL*_2.fastq.gz > $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_2.fastq.gz"
  	else # If only single file, move, don't copy
	  	mv $SEARCH_DL*_1.fastq.gz $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_1.fastq.gz"
	  	mv $SEARCH_DL*_2.fastq.gz $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_2.fastq.gz"
	  fi
 
	else #Single-End Reads
    echo "Concatenating .fastq files..." 
    COUNT=`ls -1 $SEARCH_DL*.fastq.gz | wc -l`
    if [[ $COUNT > 1 ]] ; then # If only single file, move, don't copy
      cat $SEARCH_DL*.fastq.gz > $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME".fastq.gz"
    else
      mv $SEARCH_DL*.fastq.gz $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME".fastq.gz"
    fi
	fi

  rm $SEARCH_DL
}

### Trim fastq files for adaptors using trimmomatic
### Could be called with $1=directory that holds fastq files
function trimReads () {
  if $TRIM_READ || [[ $1 != "" ]]; then #if a directory is fed into this function, then function will run on the directory
    cd $TEMP_DIR

    for FILE in $CURRENT_DIRECTORY/${1:-$FASTQ_DIRECTORY}/*$SEARCH_KEY*fastq.gz; do
      determinePairedFastq

      if [[ $PAIRED_END ]]; then
        echo "Trimming "$NAME"_1.fastq.gz and "$NAME"_2.fastq.gz"

        $JAVA $TRIMMOMATIC PE -threads 10 $FASTQ_PATH"_1.fastq.gz" $FASTQ_PATH"_2.fastq.gz" \
        $NAME"_trim_1.fastq.gz" $NAME"_unpaired_trim_1.fastq.gz" $NAME"_trim_2.fastq.gz" $NAME"_unpaired_trim_2.fastq.gz" \
        ILLUMINACLIP:$ILLUMINA_ADAPATORS_ALL":2:30:10" \
        SLIDINGWINDOW:4:20 \
        MINLEN:36        
        
        echo "Moving and renaming..."
        mv $NAME"_trim_1.fastq.gz" $FASTQ_PATH"_1.fastq.gz"
        mv $NAME"_trim_2.fastq.gz" $FASTQ_PATH"_2.fastq.gz"

      else #Single-End
        echo "Trimming "$NAME".fastq.gz"

        $JAVA $TRIMMOMATIC SE -threads 10 $FASTQ_PATH".fastq.gz" $FASTQ_PATH"_trim.fastq.gz" \
        ILLUMINACLIP:$ILLUMINA_ADAPATORS_ALL":2:30:10" \
        SLIDINGWINDOW:4:20 \
        MINLEN:36

		    mv $NAME"_trim.fastq.gz" $FASTQ_PATH".fastq.gz"
      fi
    done
    
    echo "Removing unpaired trim outputs..."
    rm *unpaired*
  
  else  
    echo "Reads not trimmed"

  fi
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
  
  if [[ $NAME == *Rep* ]] ; then
    X=${NAME%_*} #removing the "_Rep"
    FOLDER_NAME=${X##*_} #Removing everything before the last _ (leaving grouping identifier)
  else
    FOLDER_NAME=${NAME##*_}
  fi
}

### Align fastq files to using data-specific aligner
### $1 can be a specific dataset prefix, will only align according to that tag
function masterAlign () {
  if $BAM_INPUT; then #if input is BAM alignment is not needed - exit function
    return 0
  fi

  cd $TEMP_DIR
  
  STAR_ARGUMENTS="--genomeDir $STAR_GENOME_DIR --runThreadN $RUN_THREAD --sjdbOverhang 70 --outFilterType BySJout --twopassMode Basic --twopass1readsN 1000000000 --outSAMunmapped Within --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --readFilesCommand zcat "

  BISMARK_ARGUMENTS="--chunkmbs $BISMARK_MEM --multicore $BOWTIE_THREAD --genome $BISMARK_GENOME_DIR "
  SEARCH_KEY=${1:-$SEARCH_KEY}

  for FILE in $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/*$SEARCH_KEY*fastq.gz; do
    determinePairedFastq

    if [[ $FILE == *"RNA"* ]]; then
      alignSTAR
    elif [[ $FILE == *"RRBS"* ]] || [[ $FILE == *"BSSeq"* ]] || [[ $FILE == *"PBAT"* ]]; then
      alignBismark
    else
      if $ALLELE_SPECIFIC || $USE_BOWTIE; then
        alignBowtie2
      else
        alignBWA
      fi
    fi
    
    readyBam
  done

  cd $CURRENT_DIRECTORY
}

### Alignment of RNASeq data using STAR
function alignSTAR () {
  FILE_STAR_OUTPUT=$NAME"Aligned.out.bam"

  if $PAIRED_END; then
    echo "Aligning "$NAME"_1.fastq.gz and "$NAME"_2.fastq.gz to genome..."
    $STAR $STAR_ARGUMENTS --readFilesIn $FILE_FASTQ1 $FILE_FASTQ2 --outFileNamePrefix $NAME

  else #Single-End
    echo "Aligning $NAME.fastq.gz to genome..."
    $STAR $STAR_ARGUMENTS --readFilesIn $FILE_FASTQ --outFileNamePrefix $NAME
  fi

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
    echo "Aligning "$NAME"_1.fastq.gz and "$NAME"_2.fastq.gz to genome..."
    $BISMARK $BISMARK_ARGUMENTS -1 $FILE_FASTQ1 -2 $FILE_FASTQ2
    BISMARK_OUTPUT=${BISMARK_OUTPUT//_bismark_bt2.bam/_1_bismark_bt2_pe.bam}

  else #Single-End
    echo "Aligning $NAME.fastq.gz to genome..."
    $BISMARK $BISMARK_ARGUMENTS $FILE_FASTQ
  fi
   
  mv $BISMARK_OUTPUT $FILE_RAW_BAM
  cleanBismark
  rm *report.txt
}

### Alignment of ChIP data with BWA
function alignBWA () {
  FILE_SAM=$NAME".sam"
     
  if $PAIRED_END; then
    echo "Aligning "$NAME"_1.fastq.gz and "$NAME"_2.fastq.gz to genome..."
    $BWA mem -t $RUN_THREAD $GENOME_FILE $FILE_FASTQ1 $FILE_FASTQ2 > $FILE_SAM
  
  else #Single-End
    echo "Aligning "$NAME" to genome..."
    $BWA mem -t $RUN_THREAD $GENOME_FILE $FILE_FASTQ > $FILE_SAM
  fi

  echo "Converting to .bam file..."
  $SAMTOOLS view -bhS -@ $RUN_THREAD $FILE_SAM > $FILE_RAW_BAM
  rm $FILE_SAM
  
}
### Alignment of ChIPSeq data using BT2 (for allelic specific)
### TODO
function alignBowtie2 () {
  $BOWTIE2 -x $GENOME_DIR -p $RUN_THREAD -l $FILE_FASTQ1 $FILE_FASTQ2 --local -S $FILE"_raw.sam"
}

### Sort and mark duplicates in BAM file after alignment
function readyBam () {
  FILE_SORTED_BAM=$NAME"_sort.bam"

  echo "Sorting bam..." #sort BAM by coordinates
  $SAMTOOLS sort -@ $RUN_THREAD -m $THREAD_MEM -o $FILE_SORTED_BAM $FILE_RAW_BAM
    
  echo "Marking duplicates..."
  $JAVA $MARKDUPS I=$FILE_SORTED_BAM O=$FILE_BAM M=temp.txt
  rm temp.txt
  rm $FILE_SORTED_BAM

  mkdir -p $CURRENT_DIRECTORY/$FOLDER_NAME/"RAW_BAM"
  mv $FILE_RAW_BAM $CURRENT_DIRECTORY/$FOLDER_NAME/"RAW_BAM"/
  mv $FILE_BAM $CURRENT_DIRECTORY/$FOLDER_NAME/

}

### Searches for .bam files containing _Rep#.bam in nested directories and combines them into a single bam
function collapseReplicates () {
  cd $TEMP_DIR
	CURRENT=""
	for FILE in $CURRENT_DIRECTORY/*/*$SEARCH_KEY*.bam; do
		if [[ $FILE == *_Rep* ]] ; then 
		  MERGED_BAM=${FILE//_Rep*.bam/.bam}
			if [[ $CURRENT != ${FILE//_Rep*.bam/}*.bam ]] ; then #
				CURRENT=$FILE
				echo "Merging ${FILE//_Rep*.bam}*.bam"
				$SAMTOOLS merge -@ $RUN_THREAD $MERGED_BAM ${FILE//_Rep*.bam/}*.bam
        if [[ $KEEP_REPLICATES == false ]]; then
				  rm ${FILE//_Rep*.bam/}_*ep*.bam #removing only the replicates
        fi
        echo "Indexing BAM file..."
        $SAMTOOLS index ${MERGED_BAM//.bam/}* #index merged & replicates (if available)
			fi
    elseec
      echo "Indexing BAM file..."
      $SAMTOOLS index $FILE
		fi
	done
  cd $CURRENT_DIRECTORY
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

### Checking dependencies of the functions
function checkDependencies () {
	echo "Checking Dependencies..."
	EXIT=0
	for COMMAND in "${DEPENDENCIES[@]}"; do
		echo -e "\tChecking $COMMAND..."
		command -v $COMMAND > /dev/null 2>&1 || {
			echo -e >&2 "\t\t$COMMAND not found!"
			EXIT=1
		}
	done

	if [[ $EXIT = 1 || $DEPEND = 1 ]] ; then
		exit 1
	fi
}

### Create neccessary files for reference genome
function setGenome () {
	if [ $FASTQ_ONLY = false ] ; then
		mkdir -p Track_Hub
		printf "hub <name>\nshortLabel <short name>\nlongLabel Hub to display <fill> data at UCSC\ngenomesFile genomes.txt\nemail <email>" > ./Track_Hub/hub.txt
	fi
  
  MOUSE="Mmu"
  RAT="Rno"

  case $1 in
		"mm9") 
			mkdir -p Track_Hub/mm9
			GENOME_DIR=$GENOME_DIR/$MOUSE/$1/
			CHROM_SIZES=$GENOME_DIR/"mm9.sizes"
			GENOME_FILE=$GENOME_DIR/"mm9.fa"
			STAR_GENOME_DIR=$GENOME_DIR/"mm9-STAR"
			BISMARK_GENOME_DIR=$GENOME_DIR
			if [ $FASTQ_ONLY = false ]; then
				printf "chr5\t143666090\t143666091" > Actb.bed
				TRACK_FOLDER="./Track_Hub/mm9/"
				mkdir $TRACK_FOLDER
				TRACKDB="./Track_Hub/mm9/trackDb.txt"
				printf "genome mm9\ntrackDb mm9/trackDb.txt" > ./Track_Hub/genomes.txt
			fi
			;;
		"mm10")
		  GENOME_DIR=$GENOME_DIR/$1/ #currently BRC version
			CHROM_SIZES=$GENOME_DIR/"mm10.sizes"
			GENOME_FILE=$GENOME_DIR/"mm10.fa"
			STAR_GENOME_DIR=$GENOME_DIR/"mm10-STAR"
			BISMARK_GENOME_DIR=$GENOME_DIR
			if [ $FASTQ_ONLY = false ] ; then
				printf "chr5\t142904365\t142904366" > Actb.bed
				TRACK_FOLDER="./Track_Hub/mm10/"
				mkdir $TRACK_FOLDER
				TRACKDB="./Track_Hub/mm10/trackDb.txt"
				printf "genome mm10\ntrackDb mm10/trackDb.txt" > ./Track_Hub/genomes.txt
			fi
			;;
		"rn6")
		  GENOME_DIR=$GENOME_DIR/$RAT/$1/
			CHROM_SIZES=$GENOME_DIR/"rn6.sizes"
			GENOME_FILE=$GENOME_DIR/"rn6.fa"
			STAR_GENOME_DIR=$GENOME_DIR/"rn6-STAR"
			BISMARK_GENOME_DIR=$GENOME_DIR
			if [ $FASTQ_ONLY = false ] ; then
				printf "chr12\t13718023\t13718024" > Actb.bed
				TRACK_FOLDER="./Track_Hub/rn6/"
				mkdir $TRACK_FOLDER
				TRACKDB="./Track_Hub/rn6/trackDb.txt"
				printf "genome rn6\ntrackDb rn6/trackDb.txt" > ./Track_Hub/genomes.txt
			fi
			;;
		"rn5")
		  GENOME_DIR=$GENOME_DIR/$RAT/$1/
			CHROM_SIZES=$GENOME_DIR/"rn5.sizes"
			GENOME_FILE=$GENOME_DIR/"rn5.fa"
			STAR_GENOME_DIR=$GENOME_DIR/"rn5-STAR"
			BISMARK_GENOME_DIR=$GENOME_DIR
			if [ $FASTQ_ONLY = false ] ; then
				printf "chr12\t15748011\t15748012" > Actb.bed
				TRACK_FOLDER="./Track_Hub/rn5/"
				mkdir $TRACK_FOLDER
				TRACKDB="./Track_Hub/rn5/trackDb.txt"
				printf "genome rn5\ntrackDb rn5/trackDb.txt" > ./Track_Hub/genomes.txt
			fi
			;;
    "hg19")
      GENOME_DIR=$GENOME_DIR/"Hsa"/$1/
			CHROM_SIZES=$GENOME_DIR/"hg19.sizes"
			GENOME_FILE=$GENOME_DIR/"hg19.fa"
			STAR_GENOME_DIR=$GENOME_DIR/"hg19-STAR"
			BISMARK_GENOME_DIR=$GENOME_DIR
			if [ $FASTQ_ONLY = false ] ; then
				printf "chr7\t5527531\t5527532" > Actb.bed
				TRACK_FOLDER="./Track_Hub/hg19/"
				mkdir $TRACK_FOLDER
				TRACKDB="./Track_Hub/hg19/trackDb.txt"
				printf "genome hg19\ntrackDb hg19/trackDb.txt" > ./Track_Hub/genomes.txt
			fi
			;;
		"oryCun2")
		  GENOME_DIR=$GENOME_DIR/"oryCun"/$1/
			CHROM_SIZES=$GENOME_DIR/"oryCun2.sizes"
			GENOME_FILE=$GENOME_DIR/"oryCun2.fa"
			STAR_GENOME_DIR=$GENOME_DIR/"oryCun2-STAR"
			BISMARK_GENOME_DIR=$GENOME_DIR
			if [ $FASTQ_ONLY = false ] ; then
				printf "chr7\t87232876\t87232877" > Actb.bed
				TRACK_FOLDER="./Track_Hub/oryCun2/"
				mkdir $TRACK_FOLDER
				TRACKDB="./Track_Hub/oryCun2/trackDb.txt"
				printf "genome oryCun2\ntrackDb oryCun2/trackDb.txt" > ./Track_Hub/genomes.txt
			fi
			;;
		"mesAur1")
		  GENOME_DIR=$GENOME_DIR/"mesAur"/$1/
			CHROM_SIZES=$GENOME_DIR/"mesAur1.sizes"
			GENOME_FILE=$GENOME_DIR/"mesAur1.fa"
			STAR_GENOME_DIR=$GENOME_DIR/"mesAur1-STAR"
			BISMARK_GENOME_DIR=$GENOME_DIR
			if [ $FASTQ_ONLY = false ] ; then
				printf "KB708222.1\t2764075\t2764076" > Actb.bed
				TRACK_FOLDER="./Track_Hub/mesAur1/"
				mkdir $TRACK_FOLDER
				TRACKDB="./Track_Hub/mesAur1/trackDb.txt"
				printf "genome mesAur1\ntrackDb mesAur1/trackDb.txt" > ./Track_Hub/genomes.txt
			fi
			;;
		*)
			echo $1 "is not a valid genome build. Enter -h for help."
			exit 1
			;;
	esac
}

### Cleaning buffer files produced by STAR during RNA alignment
#function cleanSTAR () { #TODO Make this more specific?
#	rm -r $TEMP_DIR/*$SEARCH_KEY*"STAR"
#	rm $TEMP_DIR/*tab
#	rm -r $TEMP_DIR/*"STAR"*
#}

function cleanFASTQ () {
  if [[ $KEEP_FASTQ == false ]]; then
    echo "Not keeping fastq..."
    local FASTQ_DIR=$CURRENT_DIRECTORY/$FASTQ_DIRECTORY
    rm $FASTQ_DIR/*$SEARCH_KEY*fastq.gz

    if [[ $(ls -1 $FASTQ_DIR | wc -l) == 0 ]]; then #if the Fastq directory isempty, then it can be removed
      rm -r $FASTQ_DIR
    fi

  fi
}

function cleanBAMS () { # Remove noncanonical chr alignments
	TEMP=$TEMP_DIR/"temp.bam"
	for FILE in */*.bam
	do
		echo "Cleaning" $FILE
		$SAMTOOLS view -h $FILE | awk 'BEGIN {OFS="\t"} {
			if ( $1 !~ /^@/ ) {
				if ( $3 !~ /chr[0-9XYM]*$/ ) {
					$3 = "*";
					$4 = 0;
					$5 = 0;
				}
			}
			print $0;
		}' | $SAMTOOLS view -b -o $TEMP
		$SAMTOOLS sort -@ $RUN_THREAD -m $THREAD_MEM -o $FILE $TEMP
		rm $TEMP
	done
}

function fixMates () {
	TEMP=$TEMP_DIR/"temp.bam"
	SORT=$TEMP_DIR/"sort.bam"
	for FILE in */*.bam
	do
		echo "Fixing mates in" $FILE
		$SAMTOOLS sort -n -@ $RUN_THREAD -m $THREAD_MEM -o $SORT $FILE
		$SAMTOOLS fixmate $SORT $TEMP
		$SAMTOOLS sort -@ $RUN_THREAD -m $THREAD_MEM -o $FILE $TEMP
		rm $SORT
	done
}

############### TRACK-HUB FUNCTIONS ###############
function masterTrackHub () {
  BAM_COVERAGE_ARGUMENTS="--binSize $BIN_SIZE -p $RUN_THREAD --normalizeUsing $NORMALIZE --smoothLength $SMOOTH --outFileFormat bigwig --minMappingQuality $MIN_MAPPING_QUAL --ignoreDuplicates "

  PRINTED_DIR=""
  
  for FILE_1 in ./*/*${SEARCH_KEY//_Rep*/}*.bam; do
    FILE_2=${FILE_1//.\//} #getting rid of the "./"
    FILE=$(basename $FILE_1) #leaving just the basename of the file
    FILE_NAME=${FILE//.bam/_}$NORMALIZE ##basename without file extension

    if [[ $BIN_SIZE != 1 ]] ; then
		  FILE_NAME+="_b"$BIN_SIZE
	  fi
	  if [[ $SMOOTH != 0 ]] ; then
		  FILE_NAME+="_s"$SMOOTH
	  fi

    FOLDER_NAME=${FILE_2%%\/*} #removing longest text of the matching pattern

    if [[ $FOLDER_NAME != $PRINTED_DIR ]] ; then # Create supertrack for housing associated tracks
		  printf "track %s\nsuperTrack on show\nshortLabel %s\nlongLabel %s\n\n" $FOLDER_NAME $FOLDER_NAME $FOLDER_NAME | tee -a $TRACKDB 
      # %s inteprets arguments as a string (here we are basically inserting $FOLDER_NAME into the file)
		  PRINTED_DIR=$FOLDER_NAME
	  fi

	  determineInfoBAM

    if [[ $TYPE == "RNAseq" ]] ; then
		  generateRNATrack
	  elif [[ $TYPE == "Bisulfite-Seq" ]] ; then
      generateBSTrack
	  else #ChIPseq - not stranded
		  generateBigwigsUnstranded $FILE_2 $FILE_NAME
		  printTrackHubUnstranded $FOLDER_NAME $FILE_NAME
	  fi
  done
  rm Actb.bed
}

function determineInfoBAM () {
  if [[ $FILE == *"RNA"* ]] ; then
		TYPE="RNAseq"
	elif [[ $FILE == *"RRBS"* ]] || [[ $FILE == *"BSSeq"* ]] || [[ $FILE == *"PBAT"* ]] ; then
		TYPE="Bisulfite-Seq"
	else
		TYPE="ChIPseq"
	fi
  echo "Data are" $TYPE

	FLAG=$($SAMTOOLS view $FILE_2 | head -n 1 | cut -f2)
	if [[ $((FLAG&1)) == 1 ]] ; then #Paired-end
		PAIRED=true
	else
		PAIRED=false
	fi
	echo "Data are Paired-end:" $PAIRED
}

function generateRNATrack () {
  if [ $PAIRED = true ] ; then
			  echo "Extracting F reads over Actb..."
			  $SAMTOOLS view -L Actb.bed -f 64 $FILE_2 > Actb.sam
		  else
			  echo "Extracting reads over Actb..."
			  $SAMTOOLS view -L Actb.bed $FILE_2 > Actb.sam
		  fi
		  STRANDED=`awk 'BEGIN{PLUS=0; MINUS=0} {
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
		  }' Actb.sam`
		  rm Actb.sam
		  echo "Data are" $STRANDED
		  if [[ $STRANDED == "Unstranded" ]] ; then
			  generateBigwigsUnstranded $FILE_2 $FILE_NAME
			  printTrackHubUnstranded $FOLDER_NAME $FILE_NAME
		  else
			  FILE_BIGWIG_POS=$TRACK_FOLDER"$FILE_NAME""_pos.bw"
			  FILE_BIGWIG_NEG=$TRACK_FOLDER"$FILE_NAME""_neg.bw"
			  FILE_BIGWIG_TEMP=$TEMP_DIR"/temp.bw"
			  FILE_BEDGRAPH_TEMP=$TEMP_DIR"/temp.bedgraph"
			  FILE_BEDGRAPH_TEMP2=$TEMP_DIR"/temp2.bedgraph"
			  $BAMCOVERAGE $BAM_COVERAGE_ARGUMENTS -b $FILE_2 --filterRNAstrand forward --outFileName $FILE_BIGWIG_POS
			  $BAMCOVERAGE $BAM_COVERAGE_ARGUMENTS -b $FILE_2 --filterRNAstrand reverse --outFileName $FILE_BIGWIG_NEG
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

  #		FILE_BEDGRAPH=$TEMP_DIR/$FILE_NAME".bedGraph"
		  FILE_TEMP_1=$TEMP_DIR"/temp.bam"
		  TEMP_BEDGRAPH=${FILE_TEMP_1//.bam/.bedGraph}
  #		FILE_TEMP_2=$TEMP_DIR"/temp2.bam"
		  echo "Filtering" $FILE "for mapping quality..."
		  $SAMTOOLS view -bh -@ $RUN_THREAD -q $MIN_MAPPING_QUAL -o $FILE_TEMP_1 $FILE_2
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


function generateBigwigsUnstranded () { #$1=Name of filtered .bam file $2=Final name
	echo "Generating bigwig files..."
	$SAMTOOLS index $1
	$BAMCOVERAGE $BAM_COVERAGE_ARGUMENTS -b $1 --outFileName "$TRACK_FOLDER$2".bw
}

function printTrackHubUnstranded () { #$1=Supertrack name $2=Track name
	printf "\ttrack %s\n\tparent %s\n\tshortLabel %s\n\tlongLabel %s\n\ttype bigWig\n\tbigDataUrl %s\n\tcolor 200,50,0\n\tvisibility full\n\tmaxHeightPixels 100:60:25\n\tautoScale on\n\talwaysZero on\n\n" $2 $1 $2 $2 $2".bw" | tee -a $TRACKDB
}

function printTrackHubStranded () { #Takes in supertrack name then track name
	printf "\ttrack %s\n\tparent %s\n\tcontainer multiWig\n\tshortLabel %s\n\tlongLabel %s\n\ttype bigWig\n\tvisibility full\n\tmaxHeightPixels 100:60:25\n\tconfigurable on\n\tautoScale on\n\talwaysZero on\n\taggregate transparentOverlay\n\tshowSubtrackColorOnUi on\n\tpriority 1.0\n\n" $2 $1 $2 $2 | tee -a $TRACKDB
	printf "\t\ttrack %s\n\t\tparent %s\n\t\tshortLabel %s\n\t\tlongLabel %s\n\t\ttype bigWig\n\t\tbigDataUrl %s\n\t\tcolor 200,50,0\n\t\tautoScale on\n\n" $2"_pos" $2 $2"_pos" $2"_pos" $2"_pos.bw" | tee -a $TRACKDB
	printf "\t\ttrack %s\n\t\tparent %s\n\t\tshortLabel %s\n\t\tlongLabel %s\n\t\ttype bigWig\n\t\tbigDataUrl %s\n\t\tcolor 200,50,0\n\t\tautoScale on\n\n" $2"_neg" $2 $2"_neg" $2"_neg" $2"_neg.bw" | tee -a $TRACKDB
}

