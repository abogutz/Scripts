

STEPS TO SET UP PIPELINE ON SERVER
1. Download miniconda.
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
Do you wish the installer to initialize Miniconda3 by running conda init? [yes]

Restart terminal to allow for initialization of conda. 


2. Download all required software for the pipeline.
The ones listed are the ones required for Graham/Cedar.
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


4. Open the config file and complete all the TODOs.

5. Open the MasterDAT.sh and complete all the TODOs.

Also, edit the top of the script with the corresponding submission requirement. 
Fill in any <> spaces with your own information.

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
