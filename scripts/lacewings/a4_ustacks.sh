#!/bin/bash
#SBATCH --job-name=ustacks
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-32]%20

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

module load stacks/2.41

# this is an array job. the code will be run 33 times, once for each sample, with 20 jobs running simultaneously.
# each instance will have a SLURM_ARRAY_TASK_ID between 0 and 32. 
# we use the SLURM_ARRAY_TASK_ID to grab a different sample for each job. 

#input/output directories, supplementary files
INDIR=../../results/demultiplexed_fastqs

# make output directory if it doesn't exist
OUTDIR=../../results/stacks/denovo
mkdir -p $OUTDIR

# make an array of all fastq files
FASTQS=($(ls -1 $INDIR/*fq.gz))
# select the fastq file for this array task using SLURM_ARRAY_TASK_ID
INFILE=$(echo ${FASTQS[$SLURM_ARRAY_TASK_ID]})

# specify the unique integer (-i $ID), add 1 b/c the task ids start at 0. 
ID=$(expr 1 + $SLURM_ARRAY_TASK_ID)
# specify the sample ID (--name $SAM)
# pull the sample ID from the fastq file name using grep
SAM=$(echo $INFILE | grep -oP "[0-9]{3}[^\.]+")

# run ustacks
ustacks \
-f $INFILE \
-o $OUTDIR \
-i $ID \
--name $SAM \
-t gzfastq \
-p 6 \
-M 8 \
-m 3 \
--max-gaps 10 \
--high-cov-thres 10 

