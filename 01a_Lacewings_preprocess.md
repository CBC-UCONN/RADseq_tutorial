# Demultiplexing and quality control. #

## Introduction

This section of the tutorial deals with demultiplexing and basic quality control of the lacewing dataset. 

Steps here will use the following software packages:

- [ Stacks ](http://catchenlab.life.illinois.edu/stacks/)
- [ FastQC ](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [ MultiQC ](https://multiqc.info/)

Each major step has an associated bash script tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) job scheduler. The code can easily be modified to run interactively, or using another job scheduler. 

## Contents
  
1.    [ Motivation ](#Motivation)
2.    [ Demultiplex the sample pool ](#Demultiplexing)
3.    [ Assess sequence quality with FastQC ]()
3.    [ Quality trim using Trimmomatic ]()
3.    [ Reassess sequence quality ]()

## Motivation

Before we can run a variant-calling pipeline, we need to pre-process the data. Most RAD-seq protocols involve adding sample-specific barcodes to the DNA fragments to be sequenced and then pooling all the samples together. After sequencing, we need to separate each sample out by its barcode, a process often referred to as "demultiplexing". After that, we need to do some basic checks and quality control to evaluate our data and ensure its suitability for downstream analysis. 

HTJ generated the lacewing dataset using the "original" RADseq protocol. They digested whole genomic DNA using the restriction enzyme SbfI and then ligated adapters containing barcodes to the overhang. They then randomly sheared the fragments, size-selected them, and added an adapter to the randomly sheared end. They finally sequenced the resulting DNA fragments from the end with the restriction site. 

So what does this mean for the DNA sequences we expect to see? 

First, SbfI cuts at this motif:

<img src="/img/sbfI.png" alt="SbfI cut site" width="300"/>

##Demultiplexing

Before we start, it will be helpful to visualize 
