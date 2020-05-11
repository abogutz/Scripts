

STEPS TO SET UP PIPELINE ON SERVER
1. Download miniconda.
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
Do you wish the installer to initialize Miniconda3 by running conda init? [yes]

Restart terminal to allow for initialization of conda. 


2. Download all required software for the pipeline.
The ones listed are the ones required for Graham/Cedar.
conda activate
conda install -c bioconda entrez-direct
conda install -c bioconda deeptools
conda install -c bioconda bismark
conda install -c bioconda/label/cf201901 sra-tools
conda install -c bioconda hicup

#move into miniconda/bin for the next downloads (you could download them in a different location as well)
cd ~/miniconda/bin

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
chmod +x bedGraphToBigWig

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed
chmod +x bigWigAverageOverBed


3. Clone this repository from github onto server.
(Personally like to keep it at a location where I usually work in so I can access the script easily. Alternatively add it to your path.)
git clone https://github.com/tiffyl/Scripts.git


4. Open the config file and complete all the USER ACTION REQUIRED.

5. Open the MasterDAT.sh and complete all the USER ACTION REQUIRED.

6. Also, edit the top of the script with the corresponding submission requirement. 
Fill in any <> spaces with your own information. Examples:

----> COMPUTE CANADA (CEDAR/GRAHAM)
#SBATCH --account=<Group Name>            # required
#SBATCH --ntasks=8                        # number of MPI processes
#SBATCH --mem-per-cpu=12G                 # memory; default unit is megabytes
#SBATCH --time=01-12:00                   # time (DD-HH:MM)
#SBATCH --mail-user=<email address>
#SBATCH --mail-type=ALL

----> BRC
#$ -cwd				# execute script from current directory
#$ -pe ncpus 4			# parallel environment
#$ -l h_vmem=15G		# memory
#$ -m e
#$ -M <email address>



7. Check dependencies by running MasterDAT.sh with the -D option

8. Check available options by running MasterDAT.sh with the -h option

9. Run the pipeline. For example: 

Working in USER jra on the CEDAR server:
Starting with a list of SRA codes
~/bin/Scripts/MasterDAT.sh -i Leung2020.txt
Starting with zipped FASTQ files
sbatch ~/bin/Scripts/MasterDAT.sh -k -r -T -f ./Fastq -x C57BL6JxSPRET_14wk_Female_brain_ctx_RNA_Keown2017





Troubleshooting

Issue: Output files have names 'study*fastq.gz'
Problem: You forgot to add "_rep#" to your filename, causing the script to not detect your input files
Solution: Append "_rep#" to your filename, even if there is only one replicate available, e.g.:
myfile_RNAseq_study2020.fastq.gz --> myfile_RNAseq_study2020_rep1.fastq.gz 
Problem2: You input FASTQs that did not end in exactly ".fastq.gz".

Issue: Folder containing output files does not exist
Problem: You forgot to add "_studyAndYear" to your filename, causing the script to not create a new directory
Solution example: myfile_RNAseq_rep1.fastq.gz --> myfile_RNAseq_Batten2019_rep1.fastq.gz
