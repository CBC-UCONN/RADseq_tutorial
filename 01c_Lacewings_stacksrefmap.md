# Running Stacks in _reference mapping_ mode. #

## Introduction

This section of the tutorial deals with running `Stacks` in _reference mapping_ mode on the lacewing dataset. 

Steps here will use the following software packages:

- [ Stacks ](http://catchenlab.life.illinois.edu/stacks/)
- [ bwa ](http://bio-bwa.sourceforge.net/)
- [ samtools ](http://www.htslib.org/doc/samtools.html)

Each major step has an associated bash script tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) job scheduler. The code can easily be modified to run interactively, or using another job scheduler. 

## Contents
  
1.    [ Motivation ](#Motivation)
2.    [ Step 1:  ](#Map-to-the-reference-genome)
3.    [ Step 2:  ](#)
3.    [ Step 3: populations ](#Step-2-populations)

## Motivation

In the [previous section](/01b_Lacewings_stacksdenovo.md) we ran stacks in _de novo_ mode, meaning that we assembled our sequence data into loci, called variants, and genotyped individuals without using a reference genome. This involved several clustering steps to parse out which sequences should belong to which loci both within and across individuals. In reference mapping mode, we align all the sequences to a reference genome, and membership in a locus is implied by mapping position. A major advantage of this approach is that we don't need to specify all those clustering thresholds. Additionally, distinct loci with high sequence similarity may be easier to distinguish if they are represented in the reference genome. 

This approach has only three steps:

1. Map to a reference genome using your favorite aligner. Here we'll use `bwa mem`. 
2. Run `gstacks` to identify RAD loci and genotype individuals. 
3. Run `populations` to filter loci and reformat output for various downstream programs. 

Here we'll use the wrapper `ref_map.pl` to run both `gstacks` and `populations`. We demonstrated how to run the pipeline piecemeal in the previous section, and there is no real opportunity for parallelization (after mapping), as all individuals must be considered simultaneously when genotyping. 

## Step 1: Map to the reference genome

The first step is to align the reads for each sample to the reference genome. Here we're going to use `bwa mem`, a very commonly used short read aligner. Because each sample can be mapped independently here, we have an opportunity to parallelize this work by using an SLURM array job. As explained in the previous section, for an array job, we write one script, and SLURM will run it for each of our samples, swapping out sample names for each task. To learn more about this see our guide [here](https://github.com/CBC-UCONN/CBC_Docs/wiki/Job-arrays-on-Xanadu). 

Programs that map short reads to reference genomes generally require the genome to be _indexed_ first. We have already computed the index here, and all the index files are found in the same directory as the reference genome itself. If you need to create an index, it's fairly simple, though it may take some time if your genome is very large. 

```bash
bwa index /path/to/refgenome.fasta
```

With an indexed reference genome in hand, we'll run [this script](/scripts/lacewings/c1_bwaalign.sh). 

Several key parameters are represented by shell variables (e.g. $INFILE), which are defined at the beginning of the script (and explained [here](/00_preliminaries.md#putting-it-all-together)). The core of the script looks like this:

```bash
# run bwa mem to align, then pipe it to samtools to compress, then again to sort
bwa mem -t 4 -R $RG $REFERENCE $INDIR/$INFILE | \
samtools view -S -h -u - | \
samtools sort -T /scratch/$SAM - >$OUTDIR/$OUTFILE

# index the bam file
samtools index $OUTDIR/$OUTFILE
```

This call uses pipes (`|`) to align the reads, compress them, and sort them by genomic coordinate without writing any intermediate files. We've used `\` to break this command over several lines. 

The first line invokes the aligner, `bwa mem`. The `-t` flag indicates how many CPU threads should be used, `-R` specifies the read group, then we give the reference genome and the input fastq file. The read group flag causes each SAM record to be tagged with a read group ID (for details see the [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf)). This is not strictly necessary with `stacks`, but it is good practice, and required by other standard variant callers (`bcftools`, `gatk`, `freebayes`). 

`bwa mem` by default writes to the standard output, so we use a pipe, `|` to redirect that to `samtools view`. The `-S` flag indicates the input is in uncompressed SAM format, `-h` indicates the file header should be written, `-u` says output should be in uncompressed BAM and the final `-` indicates the input should be read from the standard input. 

The BAM formatted alignments are then piped to `samtools sort` to be coordinate sorted and then written to a file. `-T` specifies a temporary file name that can be used if it necessary during sorting (it usually is). 

Finally, we index the resulting bam files. 

We can run the script by navigating to the `scripts/lacewings/` directory and typing:

```bash
sbatch c1_bwaalign.sh
```

## Step 2: refmap.pl

Now that we have a set of bam files generated by `bwa` and `samtools`, we can run the `ref_map.pl` script. This script will execute both `gstacks` to assemble the loci and genotype individuals, and `populations` to filter and output data. We'll use [this script](/scripts/lacewings/c2_refmappl.sh). 

As usual, some key files are represented as shell variables (e.g. $POPMAP) that are defined higher in the script. The call looks like this:

```bash
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
```

The flags are as follows: `--samples` gives the directory containing the bam alignment files. `--popmap` gives the population map giving the assignment of samples to populations, as explained [here](/01b_Lacewings_stacksdenovo.md#Step-2-cstacks). `-o` gives the output directory. `-T` gives the number of CPU threads to be used. The `-X` flags are used to feed command line options to the constituent parts of the `stacks` pipeline. In this case, we don't need to alter any options in `gstacks`, so we just feed the exact same options to `populations` as we did in the [previous section](/01b_Lacewings_stacksdenovo.md#Step-6-populations). 

