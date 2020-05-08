#! /bin/bash

pushd $(dirname $0) > /dev/null
source $(pwd -P)/Graham.config
popd > /dev/null

CURRENT_DIRECTORY=$(pwd)

### Using BW files to obtain counts & average over given BED file
### $1=BED file, $2=Feature of BED file (for output suffix), $3=Output directory
function avgOverBed () {
	BEDFILE=$CURRENT_DIRECTORY/$1
	OUTPUT_SUFFIX=$2
	OUTPUT_DIR=$CURRENT_DIRECTORY/$3

	cd $TEMP_DIR
	for FILE in $CURRENT_DIRECTORY/*.bw; do
		OUT_FILE=$(basename $FILE)
		OUT_FILE=${OUT_FILE//.bw/}
		echo "Calculating $OUT_FILE average over $OUTPUT_SUFFIX..."
		$BWAVGOVERBED $FILE $BEDFILE $OUTPUT_DIR/$OUT_FILE"_AvgOver"$OUTPUT_SUFFIX".bed"
	done
	cd $CURRENT_DIRECTORY
}

### Using BAM files to obtain RPKM over bed features, list of BAM 
function bam2rpkm () {
	BED=$1
	OUTFILE=$2
	TEMP_FILE="temp_"$(basename $OUTFILE)
	
	#check if output file exist, if not - set up with bed file with of chr, start, end (can be joined with names later in R)
	if [[ -f "$OUTFILE" ]]; then
		echo $OUTFILE "exists"
	else
		awk 'BEGIN {OFS="\t"; print "chr", "start", "end";} {print $1, $2, $3}' $BED > $OUTFILE
		echo "created" $OUTFILE
	fi

	for BAM in ${@:3}; do
		FILE=$(basename $BAM)
		SAMPLE=${FILE%%.*}"_RPKM"
		echo Sample is $SAMPLE
		TOTAL_READS=$($SAMTOOLS idxstats $BAM | awk -F '\t' '{s+=$3}END{print s}')
		echo The number of reads is $TOTAL_READS

		$BEDTOOLS multicov -bams $BAM -bed $BED \
		| awk 'BEGIN {OFS="\t"; print "chr", "start", "end", "name", "'$SAMPLE'"} {print $1, $2, $3, $4, ($5*1000000000)/(($3-$2)*'$TOTAL_READS');}' \
		| cut -f5 | paste $OUTFILE - > $TEMP_FILE
		
		mv $TEMP_FILE $OUTFILE
		echo $SAMPLE "done for" $OUTFILE
	done	
}


### split BED file into 2 files based on strandedness (+ or -)
function splitStrands () {
	BED=$1
	COLUMN=$2

	awk 'BEGIN {OFS="\t"} {if ($'$COLUMN' == "+") {print $0}}' $BED > "pos_"$BED
	awk 'BEGIN {OFS="\t"} {if ($'$COLUMN' == "-") {print $0}}' $BED > "neg_"$BED
}

