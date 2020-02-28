#!/bin/bash 
#SBATCH --job-name=populations
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=15G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date

# load software
module load stacks/2.41

INDIR=../../results/stacks/denovo

POPMAP=../../metadata/popmap.txt

populations \
-P $INDIR \
-M $POPMAP \
--vcf \
--fasta-samples \
--fasta-loci \
--treemix \
-t 8