#! /bin/bash

source $FUNCTIONS_DIR/barb.config

CURRENT_DIRECTORY=$(pwd)
DIPLOID_GENOME_DIR=$GENOME_DIR/diploid #GENOME_DIR determined from the start of the script using GENOME_BUILD
HAPLOID_GENOME_DIR=$GENOME_DIR/haploid

function checkPseudogenome() {
  if $ALLELE_SPECIFIC; then

    cd $TEMP_DIR
    
    CROSS_LIST=$(awk '{if ($3!="" || $4!="") {print $3"_"$4}}' $INPUT_FILE | sort | uniq) #create array of crosses
    EXIT_SCRIPT=false

    for CROSS in $CROSS_LIST; do
      STRAIN_1=${CROSS//_*}
      STRAIN_2=${CROSS##*_}

      if [[ -z $STRAIN_1 || -z $STRAIN_2 ]]; then
        echo -e "Error:\t One of your samples only have 1 strain entered. \nPlease ensure both strains are entered for allele specific pipeline"
        EXIT_SCRIPT=true
        continue #continue to next iteration of cross
      fi

      MAKE_GENOME=true
      
      for DIP in $DIPLOID_GENOME_DIR/*; do #searching for folder with name of both strains
        if [[ $DIP == *$STRAIN_1* && $DIP == *$STRAIN_2* ]]; then
          MAKE_GENOME=false

          if [[ ! -f $HAPLOID_GENOME_DIR/$STRAIN_1/$STRAIN_1".fa.refmap" || \
                ! -f $HAPLOID_GENOME_DIR/$STRAIN_2/$STRAIN_2".fa.refmap"]; then #check RefMaps exist
            echo -e "Error:\t One of the refmap ($STRAIN_1 or $STRAIN_2) doesn't exists."
            EXIT_SCRIPT=true
          fi

          echo "Reference genome for $CROSS is found"
          break #exit the loop of iterating thru all the diploid genomes for this combo of cross 
        fi
      done

      if $MAKE_GENOME; then
        echo "Error:\t Diploid genome was not found for $CROSS. \nPlease use provided CreateDipPseudoGenome.sh to obtain diploid pseudogenome."
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
  STRAIN_1=$(grep -e $NAME $INPUT_FILE | cut -f3 | head -n 1)
  STRAIN_2=$(grep -e $NAME $INPUT_FILE | cut -f4 | head -n 1)

  CROSS_COMBO=$STRAIN_1"_"$STRAIN_2

  if [[ -z $STRAIN_1 || -z $STRAIN_2 ]]; then #if either strain is not entered, no allele specific for this data
    echo "Strains not provided, allele specific pipeline is not ran on" $NAME
    continue #dont continue the function for this SRACODE entry (move on to next code)
  fi
  
  for DIP in $DIPLOID_GENOME_DIR/*; do #searching for folder with name of both strains
    if [[ $DIP == *$STRAIN_1* && $DIP == *$STRAIN_2* ]]; then
      GENOME_FILE=$DIP/*".fa"
      STAR_GENOME_DIR=$DIP/*"-STAR"
      BOWTIE2_GENOME_DIR=$DIP

      STRAIN_1_REF=$HAPLOID_GENOME_DIR/$STRAIN_1/$STRAIN_1".fa.refmap"    
      STRAIN_2_REF=$HAPLOID_GENOME_DIR/$STRAIN_2/$STRAIN_2".fa.refmap"

      echo "Reference directory for $NAME is set to $DIP"
      break
    fi
  done
}


function unpackAllelic () { #working on bam that has already aligned to the pseudogenome (need to run this twice, once with each strain
  cd $TEMP_DIR
  local STRAIN=$1
  
  for SAM in $TEMP_DIR/$SEARCH_KEY*".sam" #should be in replicates
  
  #take SAM file, align to pseudogenome
  #split header into two, to separate the reads into two haplotypes, rename chr from hap1_chr to chr
    FILE=$SEARCH_KEY"_"$STRAIN"_"${SAM##*_} #adding Strain before _[Rr]ep
    FILE=${FILE//.sam/}
        
    $SAMTOOLS view -H $SAM \ 
    | awk '($0 ~ "'$STRAIN'") {print $0}' \ #take lines that include the strain name (strain_specific chrom & commands)
    | sed 's/'$STRAIN'_chr/chr/g' > $FILE"_raw.sam" 

    # get UNIQUELY ALIGNED READS, separate into two files, only keep reads where their mate also maps to the same chromo of the same haplotype
    # uniquely aligned defined by MIN_MAPQ
    $SAMTOOLS view $SAM -q $MIN_MAPQ \
    | awk '(($3 ~ "'$STRAIN'")&&($7 ~ "'$STRAIN'" || $7 == "*" || $7 == "=")) {print $0}' \
    | sed 's/'$STRAIN'_chr/chr/g' >> $FILE"_raw.sam"
    
    # Convert SAM to BAM
    $CLEANSAM I=$FILE"_raw.sam" O=$FILE"_cleaned.sam"
    $SAMTOOLS view -bhS -@ $RUN_THREAD $FILE"_cleaned.sam" > $FILE"_raw.bam"
    $SAMTOOLS sort -@ $RUN_THREAD -m $THREAD_MEM -o $FILE"_sorted.bam" $FILE"_raw.bam"
    $JAVA $MARKDUPS I=$FILE"_sorted.bam" O=$FILE".bam" M=$FILE"_markDupeMetrics.txt"
#    cat $FILE"_markDupeMetrics.txt" >> $SEARCH_KEY"_alignLog.txt"

    mkdir -p $CURRENT_DIRECTORY/$CROSS_COMBO #CROSS_COMBO is STRAIN_1_STRAIN_2
    mv $FILE".bam" $CURRENT_DIRECTORY/$CROSS_COMBO

    rm $SAM
    rm $FILE* #removing all the buffer sam & bam files

  done

  
#    $SAMTOOLS index $FILE".bam"
    
#    echo "$FILE flagstats" >> $SEARCH_KEY"_alignLog.txt"
#    $SAMTOOLS flagstat $FILE".bam" >> $SEARCH_KEY"_alignLog.txt"
#    printProgress "[unpacking allelic "$STRAIN" reads] finished successfully"
        

    #remove PCR duplicate aligned reads
#    printProgress "[bamToProjectedWig] started with parameters F="$FLAG
#    $SAMTOOLS view -bh -F $FLAG $FILE".bam" > $FILE"_F"$FLAG".bam"
    # use bedtools to convert bam to bedGraph
}





function checkFileExists () {
  if [[ ! -f $1]; then
    echo "Error: File $1 doesn't exists"
    exit 1
  fi
}
