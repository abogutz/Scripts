#! /bin/bash
#SBATCH --account=def-mlorincz            # required (format def-name)
#SBATCH --cpus-per-task=8                        # number of cpus
#SBATCH --mem-per-cpu=4G                 # memory; default unit is megabytes
#SBATCH --time=00-02:00                   # time (DD-HH:MM)
#SBATCH --mail-user=aaron.bogutz@ubc.ca
#SBATCH --mail-type=ALL

GROUP=$SBATCH_ACCOUNT
RUN_THREAD=$SLURM_CPUS_PER_TASK
FASTQ_DIRECTORY="Fastq"
CURRENT_DIRECTORY=$(pwd)
TEMP_DIR=$SCRATCH"/"$SLURM_JOB_ID"/"
mkdir $TEMP_DIR

mkdir $FASTQ_DIRECTORY

module load gcc/9.3.0
module load sra-toolkit/2.10.8
FASTERQDUMP="fasterq-dump"

for FILE in SRA/*/*/*.sra
do
	SRA=$(basename $FILE)
	FOLDER=${FILE%\/*}
	NAME=${FOLDER%\/*}
	NAME=${NAME##*\/}
	
	echo $FOLDER
	cd $FOLDER
	echo "Dumping fastq files..."
	$FASTERQDUMP -e $RUN_THREAD -o "$TEMP_DIR/$NAME.fastq" $SRA

	echo "Compressing fastq files..."
	for DL_FASTQ in $TEMP_DIR/*.fastq; do 
		pigz -p $RUN_THREAD $DL_FASTQ
#		COMPRESS_PID_ARRAY[$COMPRESS_COUNTER]=$!
#		((COMPRESS_COUNTER++))
	done

#	for gzipid in ${COMPRESS_PID_ARRAY[*]}; do
#		wait $gzipid #wait for all gzip to finish
#	done	
	
	echo "Moving fastq.gz files..."
	cd $TEMP_DIR
	if `ls *_1.fastq.gz 1> /dev/null 2>&1`; then
		echo "Paired..."
#		exit
		COUNT=$(ls -1 *_1.fastq.gz | wc -l) # TODO DONE? make this tolerant of multiple SRA files for a single dataset
		if [[ -f $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_1.fastq.gz" ]]; then #Check to see if a file has already been dumped
			cat *_1.fastq.gz $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_1.fastq.gz" > temp1
			chgrp $GROUP temp1
			mv temp1 $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_1.fastq.gz"
			cat *_2.fastq.gz $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_2.fastq.gz" > temp2
			chgrp $GROUP temp2
			mv temp2 $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_2.fastq.gz"
		else
			if [[ $COUNT > 1 ]] ; then
				cat *_1.fastq.gz > $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_1.fastq.gz"
				cat *_2.fastq.gz > $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_2.fastq.gz"
			else # If only single file, move, don't copy
				chgrp $GROUP *.fastq.gz
				mv *_1.fastq.gz $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_1.fastq.gz"
				mv *_2.fastq.gz $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME"_2.fastq.gz"
			fi
		fi
	else #Single-End Reads
		echo "Single End..."
		COUNT=$(ls -1 *.fastq.gz | wc -l)
		if [[ -f $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME".fastq.gz" ]]; then #File already dumped; concat
			cat *.fastq.gz $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME".fastq.gz" temp
			chgrp $GROUP temp
			mv temp $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME".fastq.gz"
		else
			if [[ $COUNT > 1 ]] ; then # If only single file, move, don't copy
				cat *.fastq.gz > $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME".fastq.gz"
			else
				chgrp $GROUP *.fastq.gz
				mv *.fastq.gz $CURRENT_DIRECTORY/$FASTQ_DIRECTORY/$NAME".fastq.gz"
			fi
		fi
	fi

	rm *.fastq.gz
	cd $CURRENT_DIRECTORY
	#rm -r $NAME
done


rm -r $TEMP_DIR
