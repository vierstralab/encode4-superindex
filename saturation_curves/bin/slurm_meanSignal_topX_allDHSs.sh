#!/bin/bash
#SBATCH --mem=30G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-4501%800
#SBATCH --partition=encode4,bigmem,pool,queue0,queue2
#SBATCH -o /net/seq/data2/projects/ENCODE4Plus/figures/adding_additional_datasets/4501_Index/slurm_logs_4501/slurm.%N.%j.%a.names.out # STDOUT
#SBATCH -e /net/seq/data2/projects/ENCODE4Plus/figures/adding_additional_datasets/4501_Index/slurm_logs_4501/slurm.%N.%j.%a.names.err # STDERR

module unload R
module load R/4.0.5

#Go into working directory
wdir=/net/seq/data2/projects/ENCODE4Plus/figures/adding_additional_datasets/4501_Index
cd ${wdir}

masterlist_dir=/net/seq/data2/projects/ENCODE4Plus/indexes/index_altius_23-09-12/output

###Make Mean Signal  File
if [ -f meanSignal.txt ]
then
	echo "File already created"
else
	cut -f5,6 ${masterlist_dir}/masterlist.only_autosomes.filtered.bed | awk '{ result = $1 / $2; printf "%.4f\n", result }' -  > meanSignal.txt
fi

###Run the subsampling R script
Rscript /home/nasi4/proj/encode4-plus/encode4-superindex/saturation_curves/code_num_new_DHSs_meanSignal_topXPerc.R "k=${SLURM_ARRAY_TASK_ID}" 



