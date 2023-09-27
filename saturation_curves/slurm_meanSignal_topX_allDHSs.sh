#!/bin/bash
#SBATCH --mem=30G
#SBATCH --cpus-per-task=1
#SBATCH --array=2-4501%500
#SBATCH --partition=encode4,bigmem,pool,queue0,queue2
#SBATCH -o slurm_logs_4501/slurm.%N.%j.%a.names.out # STDOUT
#SBATCH -e slurm_logs_4501/slurm.%N.%j.%a.names.err # STDERR

module unload R
module load R/4.0.5

###Print the Number of Individuals in that system
masterlist_dir=/net/seq/data2/projects/ENCODE4Plus/indexes/index_altius_23-07-24/filtered_masterlist

###Make Mean Signal  File
if [ -f ${masterlist_dir}/meanSignal.txt ]
then
	echo "File already created"
else
	cut -f12 ${masterlist_dir}/masterlist_DHSs_dnase-07-24_all_chunkIDs.5.blacklistfiltered.bed> ${masterlist_dir}/meanSignal.txt
fi

###Run the subsampling R script
Rscript ./code_num_new_DHSs_meanSignal_topXPerc.R "k=${SLURM_ARRAY_TASK_ID}" 



