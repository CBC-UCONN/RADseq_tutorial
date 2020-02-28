#!/bin/bash 
#SBATCH --job-name=multiqc
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --qos=general
#SBATCH --partition=general


hostname
date

# load software
module load MultiQC/1.7

multiqc ../../results/fastqc_dm

mv multiqc_report.html ../../results/fastqc_dm/
mv multiqc_data ../../results/fastqc_dm/

date