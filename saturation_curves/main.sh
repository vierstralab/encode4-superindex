#!/bin/bash

###
#Run this script from desired folder
###


export script_path=/home/nasi4/proj/encode4-plus/encode4-superindex/saturation_curves


#Check if the number of input parameters is correct
if [ $# -eq 3 ]
then
	echo "Correct number of parameters"
	echo "Make sure to run this script in desired folder"
else
	echo "Parameter1: Number of Samples in Index"
	echo "Parameter2: Either meanSignal or allDHSs"
	echo "Parameter3: Quantile Percentage"
	echo "Make sure to run this script in desired folder"
	exit 1

fi

#Read Parameters
num_samples=$1
script_option=$2
percentile=$3

#Make output/error directories
mkdir -p outdir
mkdir -p errdir


# Get the directory containing the script
SCRIPT_DIR=$(dirname "$0")

#Read the file paths file
masterlist=`awk -F'=' '{if($1 == "masterlist") print $2}' ${SCRIPT_DIR}/file_paths.txt` 
num_dhs=`tail -n +2 ${masterlist} | wc -l`
echo "Number of DHSs in Masterlist: ${num_dhs}"



#Check which R script to run
if [ ${script_option} == "allDHSs" ]
then
	echo "allDHSs"
else
	echo "meanSignal"
	
	#Extract Mean Signal from Masterlist
	tail -n +2 ${masterlist} \
	| cut -f5,6  \
	|  awk '{ result = $1 / $2; printf "%.5f\n", result }' - \
	> meanSignal.txt

	if [ -f meanSignal.txt ] && [ -s meanSignal.txt ]
	then
		sbatch --array=1-${num_samples}%3 --output=outdir/slurm.%N.%j.%a.out --error=errdir/slurm.%N.%j.%a.err ${SCRIPT_DIR}/bin/slurm_meanSignal_topX_allDHSs.sh ${SCRIPT_DIR}
	else
		echo "Need meanSignal.txt file"
	fi
fi


