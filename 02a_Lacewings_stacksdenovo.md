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

`Stacks` implements the entire _de novo_ pipeline in a single perl script `denovo_map.pl`, but here we will run each constituent part separately to demonstrate it. 

Before we get started, it's worth explaining a little bit about how `Stacks` in _de novo_ mode works. The first step is `ustacks`. This step analyzes each sample separately. Reads are first sorted into groups with identical sequences. These groups fall into two classes based on a user-defined frequency threshold: "stacks", which it is hoped represent true allelic sequences, and "secondary reads", that may simply be alleles with sequencing errors. The stacks are then assembled into "loci" if their sequence divergence falls below a user-specified threshold. The secondary reads are then mapped back to the loci, again contingent on a user-specified divergence threshold. 

With a set of putative genetic loci identified **within** individuals `cstacks` is then run in order to merge loci across individuals and create a catalog of haplotypes for each locus. 

## Step 1: ustacks

The first step in the `Stacks` is to identify loci within each sample. This means

executed separately on each sample. 



secondary reads: -N

## Step 2: cstacks


## Step 3: sstacks


## Step 4: tsv2bam


## Step 5: gstacks


## Step 6: populations


