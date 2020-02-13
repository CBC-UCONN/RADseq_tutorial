#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-32]%20

hostname
date

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

# load software
module load fastqc

#input/output directories, supplementary files

INDIR=../../results/demultiplexed_fastqs
FASTQS=($(ls -1 $INDIR/*fq.gz))

OUTDIR=../../results/fastqc_dm
mkdir -p $OUTDIR

INFILE=${FASTQS[$SLURM_ARRAY_TASK_ID]}

# run fastqc. "*fq" tells it to run on all fastq files in directory "../rawdata/"
fastqc -t 2 -o $OUTDIR $INFILE

date