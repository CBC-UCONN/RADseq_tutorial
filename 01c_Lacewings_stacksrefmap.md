# Running Stacks in _reference mapping_ mode. #

## Introduction

This section of the tutorial deals with running `Stacks` in _reference mapping_ mode on the lacewing dataset. 

Steps here will use the following software packages:

- [ Stacks ](http://catchenlab.life.illinois.edu/stacks/)
- [ Stacks ](http://catchenlab.life.illinois.edu/stacks/)

Each major step has an associated bash script tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) job scheduler. The code can easily be modified to run interactively, or using another job scheduler. 

## Contents
  
1.    [ Motivation ](#Motivation)
2.    [ Step 1:  ](#Map-to-the-reference-genome)
3.    [ Step 2:  ](#)
3.    [ Step 3: populations ](#Step-2-populations)

## Motivation

In the [previous section](/01b_Lacewings_stacksdenovo.md) we ran stacks in _de novo_ mode, meaning that we assembled our sequence data into loci, called variants, and genotyped individuals without using a reference genome. This involved several clustering steps to parse out which sequences should belong to which loci both within and across individuals. In reference mapping mode, we align all the sequences to a reference genome, then membership to a locus is implied by mapping position in the reference. A major advantage of this approach is that we don't need to specify all those clustering thresholds. Additionally, distinct loci with high sequence similarity may be easier to distinguish if they are represented in the reference genome. 

## Step 1: Map to the reference genome

## Step 2: refmap.pl


