#! /bin/bash
#SBATCH --account=<Group Name>            # required (format def-name)
#SBATCH --cpus-per-task=8                        # number of cpus
#SBATCH --mem-per-cpu=4G                 # memory; default unit is megabytes
#SBATCH --time=01-12:00                   # time (DD-HH:MM)
#SBATCH --mail-user=<email address>
#SBATCH --mail-type=ALL



# Modify these to whatever you want. Should match up with SBATCH submission parameters above.
SRA_FILE=""
GENOME="mm10"
THREADS=8
MEM_THREAD="4G"

# Many other options exist for MasterDAT.sh, run with -h option to see all of them

$HOME/projects/def-mlorincz/scripts/MasterDAT/MasterDAT.sh -i $SRA_FILE -g $GENOME -t $THREADS -m $MEM_THREAD 

#$HOME/projects/def-mlorincz/scripts/MasterDAT/MasterDAT.sh -t $THREADS -m $MEM_THREAD -f ./Fastq -g $GENOME # Fastq input - will use .fastq.gz files located in ./Fastq (can change). Name format is very important

#$HOME/projects/def-mlorincz/MasterDAT/MasterDAT.sh -t $THREADS -m $MEM_THREAD -b -g $GENOME # BAM file input - will find ALL bam files in subdirectories of current location (one level down)

#$HOME/projects/def-mlorincz/MasterDAT/MasterDAT.sh -t $THREADS -m $MEM_THREAD -i $SRA_FILE -g $GENOME -F # Only produce fastq files after downloading