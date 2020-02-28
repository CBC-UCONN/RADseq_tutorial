#!/bin/bash 
#SBATCH --job-name=cstacks
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=15G
#SBATCH --qos=general
#SBATCH --partition=general

hostname
date

# load software
module load stacks/2.41

# input, output files, directories
INDIR=../../results/stacks/denovo

cstacks \
-P $INDIR \
-M ../../metadata/lacewing_popmap.txt \
-p 10 \
--max-gaps 10 \
-n 15
