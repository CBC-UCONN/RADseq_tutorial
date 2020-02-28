#!/bin/bash
#SBATCH --job-name=process_radtags
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general

# test run on general partition 2/12/20 took 30m wall clock time. 

hostname
date

module load stacks/2.41

# input/output directories, supplementary files
# data for this tutorial stored here:
INDIR=/UCHC/PublicShare/CBC_Tutorials/RADseq_tutorial_data/

POOL=$INDIR/pool.fq.gz
BARCODES=../../metadata/henrysamples.txt

# make demultiplexed directory if it doesn't exist
OUTDIR=../../results/demultiplexed_fastqs
mkdir -p $OUTDIR

process_radtags \
-f $POOL \
-o $OUTDIR \
-b $BARCODES \
-i gzfastq \
-y gzfastq \
-e sbfI \
-c \
-q \
-t 145 \
-s 20 \
--adapter_1 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
--adapter_2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
--adapter_mm 0

date