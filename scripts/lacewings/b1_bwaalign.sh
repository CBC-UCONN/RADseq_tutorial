#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=noah.reid@uconn.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-32]%20


hostname
date

# load software
module load bwa/0.7.17
module load samtools/1.9

# input, output files and directories
INDIR=../../results/demultiplexed_fastqs/

OUTDIR=../../results/aligned
mkdir -p $OUTDIR

# indexed reference genome
REFERENCE=/UCHC/PublicShare/CBC_Tutorials/RADseq_tutorial_data/lacewing_genome/redundans_metaquast_filtered.nomt.masked.fasta

# make a bash array of fastq files
FASTQS=($(ls -1 $INDIR/*.fq.gz))

# pull out a single fastq file
INFILE=$(echo ${FASTQS[$SLURM_ARRAY_TASK_ID]} | sed 's/.*\///')
# use sed to create an output file name
OUTFILE=$(echo $INFILE | sed 's/fq.gz/bam/')

# get the sample ID and use it to specify the read group. 
SAM=$(echo $OUTFILE | sed 's/\..*//')
RG=$(echo \@RG\\tID:$SAM\\tSM:$SAM)

# run bwa mem to align, then pipe it to samtools to compress, then again to sort
bwa mem -t 4 -R $RG $REFERENCE $INDIR/$INFILE | \
samtools view -S -h -u - | \
samtools sort -T /scratch/$SAM - >$OUTDIR/$OUTFILE

date
