#!/bin/bash 
#SBATCH --job-name=refmap.pl
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
INDIR=../../results/aligned

OUTDIR=../../results/stacks/refmap
mkdir -p $OUTDIR

# popmap file
POPMAP=../../metadata/lacewing_popmap.txt

# refmap.pl -s option is broken. 
ref_map.pl \
--samples $INDIR \
--popmap $POPMAP \
-o $OUTDIR \
-T 10 \
-X "populations:-p 5" \
-X "populations:-r 2" \
-X "populations:--hwe" \
-X "populations:--vcf" \
-X "populations:--treemix" \
-X "populations:--fasta-samples" \
-X "populations:--fasta-loci"

