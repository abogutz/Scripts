#! /bin/bash

source $FUNCTIONS_DIR/barb.config

DIPLOID_GENOME_DIR=$GENOME_DIR/diploid #GENOME_DIR determined from the start of the script using GENOME_BUILD
HAPLOID_GENOME_DIR=$GENOME_DIR/haploid

function checkPseudogenome() {
  if [[ $ALLELE_SPECIFIC && $PARALLEL ]]; then

    cd $TEMP_DIR
    
    CROSS_LIST=$(awk '($3!="" || $4!="") {print $3"_"$4}' $INPUT_FILE | sort | uniq) #create array of crosses
    EXIT_SCRIPT=false

    for CROSS in $CROSS_LIST; do
      HAPLO_1=${CROSS//_*}
      HAPLO_2=${CROSS##*_}

      if [[ -z $HAPLO_1 || -z $HAPLO_2 ]]; then
        echo -e "ERROR:\t One of your samples only have 1 haplotype entered. \nPlease ensure both haplotypes are entered for allele specific pipeline" \
        | tee -a 
        EXIT_SCRIPT=true
        continue #continue to next iteration of cross
      fi

      MAKE_GENOME=true
      
      for DIP in $DIPLOID_GENOME_DIR/*; do #searching for folder with name of both haplotypes
        if [[ $DIP == *$HAPLO_1* && $DIP == *$HAPLO_2* ]]; then
          MAKE_GENOME=false

          if [[ ! -f $HAPLOID_GENOME_DIR/$HAPLO_1/$HAPLO_1".fa.refmap" || \
                ! -f $HAPLOID_GENOME_DIR/$HAPLO_2/$HAPLO_2".fa.refmap"]; then #check RefMaps exist
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
      GENOME_FILE=$DIP/*".fa"
      STAR_GENOME_DIR=$DIP/*"-STAR"
      BOWTIE2_GENOME_DIR=$DIP

      HAPLO_1_REFMAP=$HAPLOID_GENOME_DIR/$HAPLO_1/$HAPLO_1".fa.refmap"    
      HAPLO_2_REFMAP=$HAPLOID_GENOME_DIR/$HAPLO_2/$HAPLO_2".fa.refmap"

      printProgress "[setPseudogenome] Reference directory for $NAME is set to $DIP"
      break
    fi
  done
}

function unpackAllelic () { #working on bam that has already aligned to the pseudogenome (need to run this twice, once with each haplotype
  cd $TEMP_DIR
  local HAPLO=$1

  printProgress "[unpackAllelic] Unpacking for $HAPLO started at [$(date)]"
  
  for TOT_RAW_BAM in $TEMP_DIR/$SEARCH_KEY*"_raw.bam" #should be in replicates (NAME1_rep1_raw.bam)
  
  #take SAM file, align to pseudogenome
  #split header into two, to separate the reads into two haplotypes, rename chr from hap1_chr to chr
    NAME_HAPLO=$SEARCH_KEY"_"$HAPLO"_q$MIN_MAPQ"${TOT_RAW_BAM//$SEARCH_KEY/} #inserting haplotype name ( -> NAME1_HAPLO_MINMAPQ_rep1_raw.bam)
    NAME=${NAME_HAPLO//_raw.bam/} #will include _rep if applicable
    FILE_RAW_BAM=$NAME"_raw.bam"

    printProgress "[unpackAlleleic $HAPLO] Obtaining haplotype-specific header..."
    $SAMTOOLS view -H $TOT_RAW_BAM \ 
    | awk '($0 ~ "'$HAPLO'") {print $0}' \ #take lines that include the haplotype name (haplotype_specific chrom & commands)
    | sed 's/'$HAPLO'_chr/chr/g' > $FILE_RAW_BAM

    # get UNIQUELY ALIGNED READS, separate into two files, only keep reads where their mate also maps to the same chromo of the same haplotype
    # uniquely aligned defined by MIN_MAPQ
    printProgress "[unpackAlleleic $HAPLO] Obtaining haplotype-specific reads with MAPQ >= $MIN_MAPQ"
    $SAMTOOLS view $TOT_RAW_BAM -q $MIN_MAPQ \
    | awk '(($3 ~ "'$HAPLO'")&&($7 ~ "'$HAPLO'" || $7 == "*" || $7 == "=")) {print $0}' \
    | sed 's/'$HAPLO'_chr/chr/g' >> $FILE_RAW_BAM

    printProgress "[unpackAlleleic $HAPLO] Finished unpacking $TOT_RAW_BAM for $HAPLO -> $FILE_RAW_BAM"
    
    refineBAM
#    cat $FILE"_markDupeMetrics.txt" >> $SEARCH_KEY"_alignLog.txt"

    rm $TOT_RAW_BAM
  done

  printProgress "[unpackAllelic] Unpacking bams for $HAPLO completed at [$(date)]"
  cd $CURRENT_DIRECTORY
  
#    $SAMTOOLS index $FILE".bam"
    
#    echo "$FILE flagstats" >> $SEARCH_KEY"_alignLog.txt"
#    $SAMTOOLS flagstat $FILE".bam" >> $SEARCH_KEY"_alignLog.txt"
#    printProgress "[unpacking allelic "$HAPLO" reads] finished successfully"
        

    #remove PCR duplicate aligned reads
#    printProgress "[bamToProjectedWig] started with parameters F="$FLAG
#    $SAMTOOLS view -bh -F $FLAG $FILE".bam" > $FILE"_F"$FLAG".bam"
    # use bedtools to convert bam to bedGraph
}

function projectAllelic () {
  if $ALLELE_RUN; then 
    cd $TEMP_DIR

    printProgress "[projectAllelic] Started at [$(date)]"
    
    for FILE_BAM in $CURRENT_DIRECTORY/$FOLDER_NAME/*".bam"
      if [[ $FILE_BAM == *$HAPLO_1* || $FILE_BAM == *$HAPLO_2* ]] #only projecting bams that were aligned allele specifically
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

          printProgress "[projectAllelic stranded RPM] Combining stranded bedgraphs for plus strand..."
          $BEDTOOLS unionbedg -i $NAME_MAP_FLAG"_first_neg_RPM.bedGraph" $NAME_MAP_FLAG"_second_pos_RPM.bedGraph" > $NAME_MAP_FLAG"_p_tmp.bedGraph"
          rm $NAME_MAP_FLAG"_first_neg_RPM.bedGraph" $NAME_MAP_FLAG"_second_pos_RPM.bedGraph"
          awk '{OFS="\t";FS="\t"} {print $1, $2, $3, $4+$5}' $NAME_MAP_FLAG"_p_tmp.bedGraph" > $NAME_MAP_FLAG"_pos_preProject.bedGraph"
          rm $NAME_MAP_FLAG"_p_tmp.bedGraph"
          prepWigAndProject $NAME_MAP_FLAG"_RPM_pos" $NAME_MAP_FLAG"_pos_preProject.bedGraph" " plus stranded RPM"
          
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
      fi
    done

    printProgress "[projectAllelic] All haplotype-specific bams have completed projection at [$(date)]"
    cd $CURRENT_DIRECTORY
  fi
}

function prepWigAndProject () {
  local FINAL_NAME=$1
  local PRE_BEDGRAPH=$2
  local PROGRESS_APPEND=$3

  printProgress "[projectAllelic$PROGRESS_APPEND] Converting $PRE_BEDGRAPH to WIG"
  local PRE_WIG=$FINAL_NAME"_preProject.wig"
  awk '
  #prepare track parameteres
  BEGIN {
         print "track type=wiggle_0 name="'$FINAL_NAME'"\ndescription="'$PRE_BEDGRAPH'""
  }
  #only process lines with 4 fields
  NF == 4 {
           print "fixedStep chrom="$1" start="$2+1" step=1 span=1" 
           for (i = 0; i < $3-$2; i++) {
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
