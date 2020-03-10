#! /bin/bash

source $FUNCTIONS_DIR/barb.config

CURRENT_DIRECTORY=$(pwd)
DIPLOID_GENOME_DIR=$PSEUDOGENOME_DIR/Diploid/$GENOME_BUILD #GENOME_BUILD was an option from starting script
HAPLOID_GENOME_DIR=$PSEUDOGENOME_DIR/Haploid/$GENOME_BUILD

function createPseudogenomes () {
  DIPLOID_REF_GENOME=""
  HAP1_REF_GENOME=""
  HAP2_REF_GENOME=""

  for code in $CODE_ARRAY; do
    SRACODE=$code
    NAME=$(grep -e $SRACODE $INPUT_FILE | cut -f2)
    STRAIN_1=$(grep -e $SRACODE $INPUT_FILE | cut -f3)
    STRAIN_2=$(grep -e $SRACODE $INPUT_FILE | cut -f4)

    if [[ $STRAIN_1 == "" || $STRAIN_2 == "" ]]; then #if either strain is not entered, no allele specific for this data
      echo "Strains not provided, allele specific pipeline is not ran on" $NAME
      continue #dont continue the function for this entry
    fi
    
    for dip_genome in $DIPLOID_GENOME_DIR/*; do
      if [[ $dip_genome == *$STRAIN_1* && $dip_genome == *$STRAIN_2* ]]; then
        DIPLOID_REF_GENOME=$dip_genome
        echo "Reference genome is $DIPLOID_REF_GENOME"
      fi
    done

    if [[ $DIPLOID_REF_GENOME == "" ]]; then
      for hap_genome in $HAPLOID_GENOME_DIR/*; do
        if [[ $hap_genome == *$STRAIN_1* ]]; then
          echo "Haploid genome for $STRAIN_1 exists"
          HAP1_REF_GENOME=$hap_genome
        elif [[ $hap_genome == *$STRAIN_2* ]]; then
          echo "Haploid genome for $STRAIN_2 exists"
          HAP2_REF_GENOME=$hap_genome  
        fi
      done
    fi

    if [[ $DIPLOID_REF_GENOME == "" && $HAP1_REF_GENOME != "" &&  $HAP2_REF_GENOME != ""]]; then
      cat $HAP1_REF_GENOME $HAP2_REF_GENOME > $DIPLOID_GENOME_DIR/$STRAIN_1"_"$STRAIN_2".fasta"
    else
      echo -e "Cannot find fasta file for diploid or haploid genome for at least one strain.\nPlease create genome fasta file for missing strain(s) before running script again.\nExiting script..."
      exit
    fi
}

function alspecSetGenome () {


$ALEA createGenome -snps-indels-separately 
}





/brcwork/lorincz_lab/jrichardalbert/development/alea.1.j.1_bis/bin/alea createGenome -snps-indels-separately \
/brcwork/lorincz_lab/jrichardalbert/alea_mm10/genomes/no_chr/C57BL6J.fasta \
/brcwork/lorincz_lab/jrichardalbert/alea_mm10/annotations_mm10/C57BL6NJ/C57BL_6NJ.mgp.v5.snps.dbSNP142.pass.vcf.gz \
/brcwork/lorincz_lab/jrichardalbert/alea_mm10/annotations_mm10/C57BL6NJ/C57BL_6NJ.mgp.v5.indels.dbSNP142.normed.pass.vcf.gz \
C57BL_6NJ C57BL6J ./insilico_C57B6NJ

