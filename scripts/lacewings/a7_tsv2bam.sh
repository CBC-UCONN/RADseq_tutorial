#!/bin/bash 
#SBATCH --job-name=tsv2bam
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=50G
#SBATCH --qos=general
#SBATCH --partition=general


hostname
date

# load software
module load stacks/2.41

# input, output files, directories

INDIR=../../results/stacks/denovo

POPMAP=../../metadata/popmap.txt

tsv2bam -P $INDIR -M $POPMAP -t 20
