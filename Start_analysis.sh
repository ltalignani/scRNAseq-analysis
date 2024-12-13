#!/bin/bash
###################configuration slurm##############################
#SBATCH -A bwambae
#SBATCH --job-name=scRNA_analysis
#SBATCH --time=01:00:00
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

# USAGE: sbatch Start_analysis.sh

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
echo -e "Aim ___________________ Pipeline for Single-Cell RNAseq experiment, from raw data to final count matrix with QC plots along the way"
echo -e "Date __________________ 2024.11.15"
echo -e "Run ___________________ sbatch Start_analysis.sh"
echo -e "Latest Modification ___ Updated for IFB-core"

# Set working directory
workdir=$(pwd)            #$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)
max_threads="30"

echo -e "Workdir is "${workdir}

###### Call snakemake pipeline ######
echo ""
echo -e "------------------------------------------------------------------------"
echo -e "########### SNAKEMAKE PIPELINE START ###########"
echo -e "------------------------------------------------------------------------"
echo ""

echo -e "Unlocking working directory:"
echo ""

snakemake --workflow-profile profile --directory ${workdir}/ --unlock

echo ""

echo -e "Dry Run:"
echo ""


snakemake --workflow-profile profile --directory ${workdir}/ --dry-run

echo ""
echo -e "Let's Run!"
echo ""

snakemake --cores 30 --workflow-profile profile --directory ${workdir}/ --local-cores 8 --use-conda --conda-frontend conda --rerun-incomplete --keep-going

###### Create usefull graphs, summary and logs ######
echo ""
echo -e "------------------------------------------------------------------------"
echo -e "########### SNAKEMAKE PIPELINE LOGS ############"
echo -e "------------------------------------------------------------------------"
echo ""

module load graphviz/2.40.1

mkdir -p ${workdir}/results/reports/ 2> /dev/null

graph_list="dag rulegraph filegraph"
extention_list="pdf png"

for graph in ${graph_list} ; do
    for extention in ${extention_list} ; do
	snakemake --workflow-profile profile --directory ${workdir}/ --${graph} | dot -T${extention} > ${workdir}/results/reports/${graph}.${extention} ;
    done
done

snakemake --workflow-profile profile --directory ${workdir} --summary > ${workdir}/results/reports/files_summary.txt

###### End managment ######
echo ""
echo -e "------------------------------------------------------------------------"
echo -e "################## SCRIPT END ###################"
echo -e "------------------------------------------------------------------------"
echo ""

find ${workdir}/results/ -type f -empty -delete                 # Remove empty file (like empty log)
find ${workdir}/results/ -type d -empty -delete                 # Remove empty directory

time_stamp_end=$(date +"%Y-%m-%d %H:%M")                        # Get date / hour ending analyzes
elapsed_time=${SECONDS}                                         # Get SECONDS counter
minutes=$((${elapsed_time}/60))                                 # / 60 = minutes
seconds=$((${elapsed_time}%60))                                 # % 60 = seconds

echo -e "End Time ______________ ${time_stamp_end}"                                       # Print analyzes ending time
echo -e "Processing Time _______ ${minutes} minutes and ${seconds} seconds elapsed"       # Print total time elapsed

echo ""
echo -e "------------------------------------------------------------------------"
echo ""
