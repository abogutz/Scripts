#! /bin/bash

# Utility to rapidly download SRA files for subsequent fastq dumping
# Should be run as nohup ./dlSRA.sh <INPUT> &
# Ideally run from data transfer node

PROJECT_DIR="/project/def-mlorincz/"
EDIRECT_PATH="$PROJECT_DIR/scripts/utilities/edirect/"
EFETCH=$EDIRECT_PATH"efetch"
ESEARCH=$EDIRECT_PATH"esearch"
TEMP_DIR="$HOME/scratch/"
PREFETCH="$PROJECT_DIR/scripts/utilities/sratoolkit.2.11.3-ubuntu64/bin/prefetch"

INPUT_FILE=$1

function downloadReads () {
	DL_COUNTER=0

	SRACODE=$1
	NAME=$(grep -e $SRACODE $INPUT_FILE | cut -f2)
	
	echo "Starting $NAME..."
	echo "Downloading $SRACODE..."
	mkdir ./SRA/$NAME
	$PREFETCH -X 100G -O ./SRA/$NAME $SRACODE &
	DL_PID_ARRAY[$DL_COUNTER]=$! #$! = the last process that was started
	((DL_COUNTER++)) #add one to the counter so next thing added to the array will be in the next position

	for pid in ${DL_PID_ARRAY[*]}
	do
		wait $pid #wait for all downloads to finish
	done

	echo "Finished $SRACODE download."
}


echo "Starting..."
while read n; do #read command will read a line of and split them into words (2 words in files)
	declare -a m=($n) #-a for an array m to be stored based on the line (SRACODE, NAME)
	CODE_ARRAY=$CODE_ARRAY" "${m[0]} #add only SRACODE to this list
	echo ${m[0]}
done < $INPUT_FILE #feeding FILE as stdin to this while read

for CODE in $CODE_ARRAY; do
	downloadReads $CODE
done



