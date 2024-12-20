#!/bin/bash
###################configuration slurm##############################
#SBATCH -A bwambae
#SBATCH --job-name=Mkreport
#SBATCH --time=00:05:00
#SBATCH -p fast
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH --mem=8GB
#SBATCH -o Cluster_logs/%x-%j-%N.out
#SBATCH -e Cluster_logs/%x-%j-%N.err
#SBATCH --mail-user=loic.talignani@ird.fr
#SBATCH --mail-type=FAIL
###################################################################

# USAGE: sbatch Make_report.sh

# Purge caches
#echo "purging caches"
#rm -rf /shared/home/ltalignani/.cache/
#rm -rf .snakemake/
#echo "Done."

# Defining cache destination
export XDG_CACHE_HOME=/shared/projects/bwambae/


###Load module
echo ""

echo -e "Load Modules:"
echo ""

module purge
module load snakemake/8.9.0 conda

# set umask to avoid locking each other out of directories
umask 002

###### About ######
echo ""
echo -e "------------------------------------------------------------------------"
echo -e "##### ABOUT #####"
echo -e "-----------------"
echo ""
echo -e "Name __________________ Single-Cell-Analysis_Workflow"
echo -e "Author ________________ Loïc Talignani"
echo -e "Affiliation ___________ UMR_MIVEGEC"
echo -e "Aim ___________________ Make the Snakemake report of the Single-Cell RNAseq analysis pipeline"
echo -e "Date __________________ 2024.11.15"
echo -e "Run ___________________ sbatch Make_report.sh"
echo -e "Latest Modification ___ Updated for IFB-core"

# Set working directory
workdir=$(pwd)            #$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)
max_threads="30"

echo -e "Workdir is "${workdir}

echo ""
echo -e "########### SNAKEMAKE REPORT ############"
echo -e "------------------------------------------------------------------------"
echo ""

snakemake --report my_pipeline_report.html

