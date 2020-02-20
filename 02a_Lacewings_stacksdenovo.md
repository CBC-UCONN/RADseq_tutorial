# Running Stacks in _de novo_ mode. #

## Introduction

This section of the tutorial deals with running `Stacks` in _de novo_ mode on the lacewing dataset. 

Steps here will use the following software packages:

- [ Stacks ](http://catchenlab.life.illinois.edu/stacks/)
- [ Stacks ](http://catchenlab.life.illinois.edu/stacks/)

Each major step has an associated bash script tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) job scheduler. The code can easily be modified to run interactively, or using another job scheduler. 

## Contents
  
1.    [ Motivation ](#Motivation)
2.    [ Step 1: ustacks ](#Step-1:-ustacks)
3.    [ Step 2: cstacks ](#Step-2:-cstacks)
4.    [ Step 3: sstacks ](#Step-3:-sstacks)
5.    [ Step 4: tsv2bam ](#Step-4:-tsv2bam)
6.    [ Step 5: gstacks ](#Step-5:-gstacks)
7.    [ Step 6: populations ](#Step-6:-populations)

## Motivation

`Stacks` implements two methods for detecting variants in sets of sequenced samples: a _de novo_ approach, in which loci are assembled by sequence similarity, and a reference-mapping approach, in which loci are identified after alignment of all sequences to a reference genome. Here we will cover the _de novo_ approach, and demonstrate the reference mapping approach in the next section. 

`Stacks` implements the entire _de novo_ pipeline in a single perl script `denovo_map.pl`, but here we will run each constituent module separately to demonstrate it. 

Before we get started, it's worth explaining a little bit about how `Stacks` works in _de novo_ mode. The first four steps collectively identify sequences within and across individuals that correspond to homologous loci. 

Step 1 is to run `ustacks`. This step analyzes **each sample separately**. Reads are first sorted into groups with identical sequences. These groups fall into two classes based on a user-defined frequency threshold: "stacks", which it is hoped represent true allelic sequences, and "secondary reads", which it is hoped are mostly alleles with sequencing errors. If a targeted locus contains a heterozygous site, then at this stage it should be represented by two stacks. The allelic stacks are then assembled into putative loci if their sequence divergence falls below a user-specified threshold. The secondary reads are then mapped back to the loci, again contingent on a user-specified divergence threshold. 

With a set of putative genetic loci identified **within** individuals we then run `cstacks` to merge loci across individuals and create a catalog of haplotypes for each locus. 

Next, we run `sstacks`, which matches loci from individuals to the catalog created by `cstacks`. 

The fourth step is to run `tsv2bam`, which essentially sorts the data locus-wise. Even with paired-end sequencing data, all steps up to this point have used only the first read. In this step, the second reads are incorporated into bam files. At this stage we have clusters of sequences across individuals that correspond to putative genomic loci. 

In step 5, we run `gstacks`. This module does three things: first, it assembles reads into contigs, which then form a kind of reference genome; second, it aligns the reads to those contigs; and third, it calls variants and genotypes individuals against the reference contigs. 

Finally, in step 6, we run `populations`, which we can use to filter and format the data, and to calculate some basic statistics. 

For further details about the implementation of `Stacks`, see Rochette et al. (2017) an Catchen et al. (2011). 

## Step 1: ustacks

`ustacks` clusters sequences and identifies loci within each sample. This means we can parallelize it by executing it separately for each sample. This is an advantage of running the separate modules on their own. 

Here we'll run `ustacks` using [this script](). Before explaining how to parallelize it, we'll first look at the `ustacks` call and its options. 

```bash
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
```

In this code, some key files are represented as shell variables (e.g. `$INFILE`) that are defined in the script. 

secondary reads: -N

## Step 2: cstacks


## Step 3: sstacks


## Step 4: tsv2bam


## Step 5: gstacks


## Step 6: populations


