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
2.    [ Exploring the sequence data ](#Exploring_the_sequence_data)
2.    [ Demultiplex the sample pool ](#Demultiplexing)
3.    [ Assess sequence quality with FastQC ]()
3.    [ Quality trim using Trimmomatic ]()
3.    [ Reassess sequence quality ]()

## Motivation

Before we can run a variant-calling pipeline, we need to pre-process the data. Most RAD-seq protocols involve adding sample-specific barcodes to the DNA fragments to be sequenced and then pooling all the samples together. After sequencing, we need to separate each sample out by its barcode, a process often referred to as "demultiplexing". After that, we need to do some basic checks and quality control to evaluate our data and ensure its suitability for downstream analysis. 

## Exploring the sequence data

Before we start, it will be helpful to recap the laboratory method, and visualize our expected sequence data. 

HTJ generated the lacewing dataset using the "original" RADseq protocol. They digested whole genomic DNA using the restriction enzyme SbfI and then ligated adapters containing 5bp barcodes to the overhang. They then randomly sheared the fragments, size-selected them, and added an adapter to the randomly sheared end. They finally sequenced the resulting DNA fragments from the end with the restriction site. 

So what does this mean for the DNA sequences we expect to see? 

First, each sequence should begin with the 5bp sample barcode sequence, for example `AACCA`. 

Next, because SbfI cuts at this palindromic motif:

<img src="/img/sbfI.png" alt="SbfI cut site" width="300"/>

we expect to see the sequence `TGCAGG`. Following that, sequences should be highly variable. 

So the first 11 bases of each sequence should look like this:

<img src="/img/seq_start.png" alt="Starting sequence" width="315"/>

Here, sample-specific 5bp barcodes represented by orange `X's`, the 6bp remainder of the SbfI cut site are red nucleotides, and the following template sequence as a purple `...` 

We can test this assumption by exploring the raw data a little bit. 

We haven't talked about the data format, yet, so first we'll do that. Our sequence data are in [fastq format](https://en.wikipedia.org/wiki/FASTQ_format). In fastq, each sequence is represented by 4 lines. The first line is the sequence name and always begins with ">". The second is the nucleotide sequence. The third is a comment line, which always begins with "+" and in almost all cases, is otherwise blank. The fourth line contains [phred-scaled base qualities](https://en.wikipedia.org/wiki/Phred_quality_score), which represent the confidence the sequencer has in its base calls. They are encoded in ASCII characters. 

The fastq files are also compressed using gzip, which is denoted by the suffix ".gz". We recommend always keeping large sequence files compressed, as nearly all relevant programs can deal with gzip-compressed files, and sequence data takes up a huge amount of storage on the Xanadu computer cluster. 

We can inspect the fastq file using `less`, a text file viewer. Navigate to the `data/` directory and enter the following:

```bash
less pool.fq.gz
```

Type `q` to exit. 

We can get a little better sense of what the sequences look like by summarizing a lot of them. From the `data/` directory enter the following:

```bash
zcat pool.fq.gz | grep -P '^[ACGTN]{5,}' | head -n 1000 | cut -c 1-11 | sort | uniq -c | sort -g
```

The result is set of counts of 11bp sequences from the first 1000 fastq records. The last ten lines should look like this:

```
      5 AACCATGCAGT
      5 AACCATGCATT
      6 AACCATGCACT
      7 AACCATGCATC
      8 AACCATGCAAG
      9 AACCATGCACG
     11 AACCATGCAGC
     62 AACCATGCACA
     63 AACCAAAGATC
    686 AACCATGCAGG
```

This result is a bit unrealistic because this is a *synthetic* pool. The sequences are sorted by sample, so these first 1000 sequences all have the barcode `AACCA`. Nevertheless, you can see that only 686 sequences contain the expected `AACCATGCAGG`. The rest of the sequences have errors in the restriction site, or are contaminant sequences that have been erroneously ligated. In a real pool, sample barcodes will be shuffled throughout the file, and you will see errors in the barcodes as well. 


