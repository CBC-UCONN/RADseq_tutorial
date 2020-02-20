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

#input/output directories, supplementary files
INDIR=../../results/demultiplexed_fastqs

FASTQS=($(ls -1 $INDIR/*fq.gz))

# make output directory if it doesn't exist
OUTDIR=../../results/stacks/denovo
mkdir -p $OUTDIR

INFILE=$(echo ${FASTQS[$SLURM_ARRAY_TASK_ID]})

# run ustacks

ID=$(expr 1 + $SLURM_ARRAY_TASK_ID)
SAM=$(echo $INFILE | grep -oP "[0-9]{3}[^\.]+")

# set min coverage high

ustacks \
-f $INFILE \
-o $OUTDIR \
-i $ID \
--name $SAM \
-M 8 \
-m 3 \
-p 6 \
--max-gaps 10 \
--high-cov-thres 10 \
-t gzfastq

