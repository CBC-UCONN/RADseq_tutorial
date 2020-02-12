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

The lacewing dataset was generated using the "original" RADseq protocol. Whole genomic DNA was digested using the restriction enzyme SbfI and then randomly sheared. Fragments were size-selected and then sequenced from the end of the fragment containing the restriction site. 
<span style="color: red;">text</span>

```
5'...<span style="color: red;">CCTGCA</span>GG...3'
3'...<span style="color: blue;">GG</span>ACGTCC...5'
```

##Demultiplexing

Before we start, it will be helpful to visualize 

