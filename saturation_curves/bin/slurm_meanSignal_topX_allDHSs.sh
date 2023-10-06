#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --partition=encode4,bigmem,pool,queue0,queue2

#Prepare R environment
module unload R
module load R/4.0.5

#Read in Parameters
script_path=$1
percentile=$2
num_samples=$3
binary_mtx_path=$4

###Run the subsampling R script
Rscript ${script_path}/bin/code_num_new_DHSs_meanSignal_topXPerc.R "k=${SLURM_ARRAY_TASK_ID}" "${percentile}" "${num_samples}" "${binary_mtx_path}"


