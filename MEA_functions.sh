#! /bin/bash

source $FUNCTIONS_DIR/barb.config

CURRENT_DIRECTORY=$(pwd)
DIPLOID_GENOME_DIR=$GENOME_DIR/diploid #GENOME_DIR determined from the start of the script using GENOME_BUILD
HAPLOID_GENOME_DIR=$GENOME_DIR/haploid

function checkPseudogenome() {
  if $ALLELE_SPECIFIC; then

    CROSS_LIST=$(awk '{if ($3!="" || $4!="") {print $3"_"$4}}' $INPUT_FILE | sort | uniq) #create array of crosses
    EXIT_SCRIPT=false

    for CROSS in $CROSS_LIST; do
      STRAIN_1=${CROSS//_*}
      STRAIN_2=${CROSS##*_}

      if [[ -z $STRAIN_1 || -z $STRAIN_2 ]]; then
        echo "One of your samples only have 1 strain entered. \nPlease ensure both strains are entered for allele specific pipeline"
        EXIT_SCRIPT=true
        continue #continue to next iteration of cross
      fi

      MAKE_GENOME=true
      
      for DIP in $DIPLOID_GENOME_DIR/*; do #searching for folder with name of both strains
        if [[ $DIP == *$STRAIN_1* && $DIP == *$STRAIN_2* ]]; then
          MAKE_GENOME=false
          checkFileExits $HAPLOID_GENOME_DIR/$STRAIN_1/$STRAIN_1".fa.refmap"
          checkFileExits $HAPLOID_GENOME_DIR/$STRAIN_2/$STRAIN_2".fa.refmap"

          echo "Reference genome for $CROSS is found"
          break #exit the loop of iterating thru all the diploid genomes for this combo of cross 
        fi
      done

      if $MAKE_GENOME; then
        echo "Diploid genome was not found for $CROSS. \nPlease use provided CreateDipPseudoGenome.sh to obtain diploid pseudogenome."
        EXIT_SCRIPT=true
      fi
      
    done

    if $EXIT_SCRIPT; then
      exit #exit script to ensure input is formatted correctly before running the whole script
    fi
  fi  
}

function setPseudogenome () {
  for code in $CODE_ARRAY; do
    SRACODE=$code
    NAME=$(grep -e $SRACODE $INPUT_FILE | cut -f2)
    STRAIN_1=$(grep -e $SRACODE $INPUT_FILE | cut -f3)
    STRAIN_2=$(grep -e $SRACODE $INPUT_FILE | cut -f4)

    if [[ -z $STRAIN_1 || -z $STRAIN_2 ]]; then #if either strain is not entered, no allele specific for this data
      echo "Strains not provided, allele specific pipeline is not ran on" $NAME
      continue #dont continue the function for this SRACODE entry (move on to next code)
    fi
    
    for DIP in $DIPLOID_GENOME_DIR/*; do #searching for folder with name of both strains
      if [[ $DIP == *$STRAIN_1* && $DIP == *$STRAIN_2* ]]; then
        GENOME_FILE=$DIP/*.fa
        STAR_GENOME_DIR=$DIP/*"-STAR"
        BOWTIE2_GENOME_DIR=$DIP

        STRAIN_1_REF=$HAPLOID_GENOME_DIR/$STRAIN_1/$STRAIN_1".fa.refmap"    
        STRAIN_2_REF=$HAPLOID_GENOME_DIR/$STRAIN_2/$STRAIN_2".fa.refmap"

        echo "Reference genome for $NAME is set to $DIP"
        break
      fi
    done

  done
}

function allespecUnpack () { #working on bam that has already aligned to the pseudogenome (need to run this twice, once with each strain
  for FILE in 
  #take SAM file, align to pseudogenome
  #split header into two, to separate the reads into two haplotypes, rename chr from hap1_chr to chr
    $SAMTOOLS view -H "$FILE".bam \
    | awk -v strain="$STRAIN" '($0 ~ strain) {print $0}' \
    | sed 's/'"$STRAIN"'_chr/chr/g' > "$FILE"_"$STRAIN"_raw.sam

    # get UNIQUELY ALIGNED READS, separate into two files, only keep reads where their mate also maps to the same chromo of the same haplotype
    # uniquely aligned defined by MIN_MAPQ
    $SAMTOOLS view "$FILE"_raw.bam \
    | awk -v strain="$STRAIN" '(($3 ~ strain)&&($7 ~ strain || $7 == "*" || $7 == "=")&&($5>='"$MIN_MAPQ"')) {print $0}' \
    | sed 's/'"$STRAIN"'_chr/chr/g' >> "$FILE"_"$STRAIN"_raw.sam
    # Convert SAM to BAM 
        $CLEANSAM I="$FILE"_"$STRAIN"_raw.sam O="$FILE"_"$STRAIN"_cleaned.sam
        rm "$FILE"_"$STRAIN"_raw.sam
        $SAMTOBAM I="$FILE"_"$STRAIN"_cleaned.sam O="$FILE"_"$STRAIN"_unsorted.bam
        rm "$FILE"_"$STRAIN"_cleaned.sam
        $SORTSAM SO=coordinate I="$FILE"_"$STRAIN"_unsorted.bam O="$FILE"_"$STRAIN"_sorted.bam
        $MARKDUPES I="$FILE"_"$STRAIN"_sorted.bam O="$FILE"_"$STRAIN".bam M="$FILE"_"$STRAIN"_markDupeMetrics.txt REMOVE_DUPLICATES=false
        cat "$FILE"_"$STRAIN"_markDupeMetrics.txt >> "$FILE"_alignLog.txt
        $SAMTOOLS index "$FILE"_"$STRAIN".bam
        echo ""$FILE"_"$STRAIN" flagstats" >> "$FILE"_alignLog.txt
        $SAMTOOLS flagstat "$FILE"_"$STRAIN".bam >> "$FILE"_alignLog.txt
        printProgress "[unpacking allelic "$STRAIN" reads] finished successfully"

        # remove PCR duplicate aligned reads
        printProgress "[bamToProjectedWig] started with parameters F="$FLAG""
        $SAMTOOLS view -bh -F "$FLAG" "$FILE"_"$STRAIN".bam > "$FILE"_"$STRAIN"_F"$FLAG".bam
        # use bedtools to convert bam to bedGraph
}

function checkFileExists () {
  if [[ ! -f "$1"]; then
    echo "Error: File $1 doesn't exists"
    exit 1
  fi
}
