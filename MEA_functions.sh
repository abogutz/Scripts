#! /bin/bash

source $FUNCTIONS_DIR/barb.config

CURRENT_DIRECTORY=$(pwd)
DIPLOID_GENOME_DIR=$GENOME_DIR/diploid #GENOME_DIR determined from the start of the script using GENOME_BUILD
HAPLOID_GENOME_DIR=$GENOME_DIR/haploid

function setPseudogenomes () {
  NEED_HAPLOID=true
  HAP1_REFGENOME=""
  HAP2_REFGENOME=""

  for code in $CODE_ARRAY; do
    SRACODE=$code
    NAME=$(grep -e $SRACODE $INPUT_FILE | cut -f2)
    STRAIN_1=$(grep -e $SRACODE $INPUT_FILE | cut -f3)
    STRAIN_2=$(grep -e $SRACODE $INPUT_FILE | cut -f4)

    if [[ -z $STRAIN_1 || -z $STRAIN_2 ]]; then #if either strain is not entered, no allele specific for this data
      echo "Strains not provided, allele specific pipeline is not ran on" $NAME
      continue #dont continue the function for this SRACODE entry (move on to next code)
    fi
    
    for dip_strains_dir in $DIPLOID_GENOME_DIR/*; do #searching for folder with name of both strains
      if [[ $dip_strains_dir == *$STRAIN_1* && $dip_strains_dir == *$STRAIN_2* ]]; then
        GENOME_FILE=$dip_strains_dir/*.fa
        STAR_GENOME_DIR=$dip_strains_dir/$STRAIN_1"_"$STRAIN_2"-STAR"
        BOWTIE2_GENOME_DIR=$dip_strains_dir
        NEED_HAPLOID=false

        checkFileExits $HAPLOID_GENOME_DIR/$STRAIN_1".fa.refmap"
        checkFileExits $HAPLOID_GENOME_DIR/$STRAIN_2".fa.refmap"

        echo "Reference genome for $NAME is found"
      fi
    done

    if $NEED_HAPLOID; then
      for hap_genome in $HAPLOID_GENOME_DIR/*.fa; do
        if [[ $hap_genome == *$STRAIN_1*]]; then
          echo "Haploid genome for $STRAIN_1 exists"
          HAP1_REFGENOME=$hap_genome
          checkFileExits $hap_genome".refmap"
          
        elif [[ $hap_genome == *$STRAIN_2* ]]; then
          echo "Haploid genome for $STRAIN_2 exists"
          HAP2_REFGENOME=$hap_genome
          checkFileExits $hap_genome".refmap"
        fi
      done
    fi

    if [[ $NEED_HAPLOID && $HAP1_REFGENOME != "" &&  $HAP2_REFGENOME != ""]]; then
      mkdir -p $DIPLOID_GENOME_DIR/$STRAIN_1"_"$STRAIN_2
      cat $HAP1_REFGENOME $HAP2_REFGENOME > $DIPLOID_GENOME_DIR/$STRAIN_1"_"$STRAIN_2/$STRAIN_1"_"$STRAIN_2".fasta"
      REFGENOME=$STRAIN_1"_"$STRAIN_2".fasta"
      echo "Reference genome is $REFGENOME"
      GENOME_FILE=$dip_strains_dir/*.fa
      STAR_GENOME_DIR=$dip_strains_dir/$STRAIN_1"_"$STRAIN_2"-STAR"
      BOWTIE2_GENOME_DIR=$dip_strains_dir
            
    else
      echo -e "Cannot find fasta file for diploid or haploid genome for at least one strain.\nPlease create genome fasta file (using alea createGenome) for missing strain(s) before running script again.\nExiting script..."
      exit
    fi

    
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
