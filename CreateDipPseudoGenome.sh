#! /bin/bash
#$ -cwd
#$ -pe ncpus 5
#$ -l h_vmem=12G
#$ -m e
#$ -M 

pushd $(dirname $0) > /dev/null
$SCRIPTS_DIR=$(pwd -P)
popd > /dev/null

source $SCRIPTS_DIR/BRC.config

#TODO: edit genome path with addition of the genome build 
GENOME_DIR=$GENOME_DIR/mm10/
DIPLOID_GENOME_DIR=$GENOME_DIR/diploid
HAPLOID_GENOME_DIR=$GENOME_DIR/haploid
STRAIN_LIST="129S1 B6NJ C57BL6J SvImJ"

OPTIONS="1:2:g:"

if ( ! getopts $OPTIONS opt); then
	echo -e "USAGE:\t $(basename $0) -1 <STRAIN_1> -2 <STRAIN_2>"
	exit 1
fi

while getopts $OPTIONS opt; do
	case $opt in
		1) 
			STRAIN_1=${OPTARG}
			;;
		2) 
			STRAIN_2=${OPTARG}
			;;
		g)
			GENOME_DIR=${OPTARG}
			;;
		\?)
			echo -e "ERROR: Invalid Option!" >&2
			exit 1
			;;
	esac
done

#Check to ensure both strains are entered
if [[ -z $STRAIN_1 || -z $STRAIN_2 ]]; then
	echo "Please enter 2 strains to create pseudogenome."
	exit
fi

#Ensure the strains are ones that we support
FAIL_STRAIN_1==true
FAIL_STRAIN_2==true

for strain in $STRAIN_LIST; do
	if [[ $FAIL_STRAIN_1 && $strain == $STRAIN_1 ]];
		FAIL_STRAIN_1=false
		echo -e "$STRAIN_1 supported."
	elif [[ $FAIL_STRAIN_2 && $strain == $STRAIN_2 ]];
		FAIL_STRAIN_2=false
		echo -e "$STRAIN_2 supported."
	fi
done

if [[ $FAIL_STRAIN_1 || $FAIL_STRAIN_2 ]]
	echo -e "At least one of your strains is not found. \nPlease double check spelling of the strains. \nIf the strain is not supported, please contact..."
	exit
fi

#Ensure we are not creating the same diploid pseudogenome
for DIP in $DIPLOID_GENOME_DIR/*; do #searching for folder with name of both strains
	if [[ $DIP == *$STRAIN_1* && $DIP == *$STRAIN_2* ]]; then
		echo "Reference diploid pseudogenome is found"
		exit 
	fi
done

#Find the corresponding haploid fasta files
for HAP in $HAPLOID_GENOME_DIR/*.fa; do
	if [[ $HAP == *$STRAIN_1*]]; then
		HAP1_REFGENOME=$HAP
	elif [[ $HAP == *$STRAIN_2* ]]; then
		HAP2_REFGENOME=$HAP
	fi
done

#Create diploid pseudogenome
DIP_NAME=$STRAIN_1"_"$STRAIN_2
DIP_REFGENOME=$DIP_NAME".fa"
mkdir -p $DIPLOID_GENOME_DIR/$DIP_NAME
cat $HAP1_REFGENOME $HAP2_REFGENOME > $DIPLOID_GENOME_DIR/$DIP_NAME/$DIP_REFGENOME

#Bowtie2 Indexing
$BOWTIE2"-build" --threads $RUN_THREAD $DIPLOID_GENOME_DIR/$DIP_NAME/$DIP_REFGENOME $DIP_NAME

#Bismark Indexing
$BISMARK"_genome_preparation" --path_to_aligner $BOWTIE2 $DIPLOID_GENOME_DIR/$DIP_NAME/

#STAR Indexing
mkdir -p $DIPLOID_GENOME_DIR/$DIP_NAME/$DIP_NAME"-STAR"
$STAR --runThreadN $RUN_THREAD --runMode genomeGenerate --genomeDir $DIPLOID_GENOME_DIR/$DIP_NAME/$DIP_NAME"-STAR" --genomeFastaFiles $DIPLOID_GENOME_DIR/$DIP_NAME/$DIP_REFGENOME

