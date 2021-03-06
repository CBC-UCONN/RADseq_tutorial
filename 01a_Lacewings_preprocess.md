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
2.    [ Exploring the sequence data ](#Exploring-the-sequence-data)
3.    [ Demultiplex the sample pool ](#Demultiplex-the-sample-pool)
4.    [ Assess sequence quality with FastQC ](Assess_sequence_quality_with_FastQC)

## Motivation

Before we can run a variant-calling pipeline, we need to pre-process the data. Most RAD-seq protocols involve adding sample-specific barcodes to the DNA fragments to be sequenced and then pooling all the samples together. After sequencing, we need to separate each sample out by its barcode, a process often referred to as "demultiplexing". After that, we need to do some basic checks to evaluate our data and ensure its suitability for downstream analysis. 

## Exploring the sequence data

Before we start, it will be helpful to recap the laboratory method, and visualize our expected sequence data. 

HTJ generated the lacewing dataset using the "original" RADseq protocol. They digested whole genomic DNA using the restriction enzyme SbfI and then ligated adapters containing 5bp barcodes to the overhang. Then they randomly sheared the fragments, size-selected them, and added an adapter to the randomly sheared end. They finally sequenced the resulting DNA fragments from the end with the restriction site. 

So what does this mean for the DNA sequences we expect to see? 

First, each sequence should begin with the 5bp sample barcode sequence, for example `AACCA`. 

Next, because SbfI cuts at this 8bp palindromic motif:

<img src="/img/sbfI.png" alt="SbfI cut site" width="300"/>

and separates the red from the gray fragment, we expect to see the sequence `TGCAGG`. Following that, sequences should be highly variable. 

So the first 11 bases of each sequence should look like this:

<img src="/img/seq_start.png" alt="Starting sequence" width="315"/>

Here, sample-specific 5bp barcodes represented by orange `X`'s, the 6bp remainder of the SbfI cut site as red nucleotides, and the following template sequence as a purple `...` 

We can test this assumption by exploring the raw data a little bit. 

We haven't talked about the data format, yet, so first we'll do that. Our sequence data are in [fastq format](https://en.wikipedia.org/wiki/FASTQ_format). In fastq, each sequence is represented by 4 lines. The first line is the sequence name and always begins with "@". The second is the nucleotide sequence. The third is a comment line, which always begins with "+" and in almost all cases, is otherwise blank. The fourth line contains [phred-scaled base qualities](https://en.wikipedia.org/wiki/Phred_quality_score), which represent the confidence the sequencer has in its base calls. They are [encoded as ASCII characters](https://drive5.com/usearch/manual/quality_score.html). 

Here are 3 example fastq records:

```
@70_8_1101_13464_1191/1
AACCATGCAGGCGCTTCACCTTCACCTTGCGCTTCTTGGCACGGAGGACATGGTCCTTGTGCTTCTCCAGCTGACCGAGGTGCTCGCGGCTCTTGGGCTGGGAGCGCTCCTGGTGCGTCTTGCGCTTCAGGTGCTTGTCCAGCCTGCTCGC
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJAJJJJJJJFJJJJJJJJJJJJJJFFFFJJJFJJJJJFJJJJJJJFJJFFJJJJ77A<JJJJJJ<JF<AJFJAAA<FJJJJJAFJFJFJJJJJJ-FF7
@70_8_1101_28280_1191/1
AACCATGCAGGCTCCTCCTTGCCCCGCTGCAGCGGGCCTGGCTGCGTGCATTGCTCTTCCTTATCGGGCACATTCTCCCCGACCAGCCACCAACTACCATTGCACAGCCGCACCGCACTGGCGGGGGTGGCACCCCGACCACCGCCGCCGT
+
JJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJAJFJJJJJJFJJJJJJJJJJJFJJJJJJJJFJJJJJJFJFFJJJJAJJJFJJFFJFF<JJJJJJJJJJJ-JJJ-<JJJJJ)<F-<)7AFJJJ<FJJ)
@70_8_1101_4838_1226/1
AACCATGCAGGCGCTGCGTCATGCTGCCACTGGCACACACGACGAGGTGCAGCGGCGGTGCCGTTGCGTCGAGGTTGACGCTCTGCCGGTAGACGCGCTCGCGGACGCCCCAAAGGAAGTAGGAGAGGAGAGGATTCATGCGGTTGAGGTT
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJFAFJJJJJJFFJJJJJJJJJJJJ-FJJFJAFJFJJFF<FJF-AJ7AAF-A-<<<AFFJF-<FJ<FF<A))7-<F7-
```

The fastq files are also compressed using gzip, which is denoted by the suffix ".gz". We recommend always keeping large sequence files compressed, as nearly all relevant programs can deal with gzip-compressed files, and sequence data takes up a huge amount of storage on the Xanadu computer cluster. 

We can inspect the fastq file using `less`, a text file viewer. Navigate to the `data/` directory and enter the following:

```bash
less lacewing_pool.fq.gz
```

Type `q` to exit. 

We can get a little better sense of what the sequences look like by summarizing a lot of them. From the `data/` directory enter the following:

```bash
zcat lacewing_pool.fq.gz | grep -P '^[ACGTN]{5,}' | head -n 1000 | cut -c 1-11 | sort | uniq -c | sort -g
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


## Demultiplex the sample pool

Now that we understand how the data are structured, we can proceed to demultiplexing the sample pool. In the directory `data/` is a file, `pool.fq.gz`, this is our sample pool. We're going to use the module `process_radtags` from `Stacks` to do the demultiplexing. 

We will use the script [a1_process_radtags.sh](/scripts/lacewings/a1_process_radtags.sh) to accomplish this. It is located in the `scripts/lacewings/` directory. The core of the script is below:

```bash
process_radtags \
-f $POOL \
-o $OUTDIR \
-b $BARCODES \
-i gzfastq \
-y gzfastq \
-e sbfI \
-c \
-q \
-t 145 \
-s 20 \
--adapter_1 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
--adapter_2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
--adapter_mm 0
```

In this code, some key files are represented as shell variables (e.g. `$POOL`) that are defined in the script. The flags are as follows: `-f` gives the sample pool. `-o` gives the output directory. `-b` is the barcode file. The file relating samples to barcodes is here `metadata/henrysamples.txt`. `-i` and `-y` indicate the input and output should be gzipped fastq files. `-e` indicates the restriction enzyme is SbfI. `-c` and `-q` remove reads with N's and overall low quality respectively. `-t` truncates reads to 145bp from the original 150bp (in these data the quality drops in the last few cycles). `-s` removes reads whose average quality in any window (0.15x read length) drops below 20. The `--adapter` flags give the adapter sequences and a mismatch tolerance, so that reads with adapter sequence can be filtered out. 

It is worth noting that if the barcode + restriction enzyme recognition sequence comes within the mismatch tolerance of the adapter sequence, you can lose **all** the reads from that sample due to the adapter filtering. This happens with one sample in this dataset if you specify `--adapter_mm 2`. 

To understand the Slurm header in the script, see our documentation in the [Xanadu tutorial](https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/xanadu/) and our [guidance on resource requests](https://github.com/CBC-UCONN/CBC_Docs/wiki/Requesting-resource-allocations-in-SLURM). 

We can run the script by navigating to the `scripts/lacewings/` directory and typing

```bash
sbatch a1_process_radtags.sh
```

This script takes 30 minutes to run. `process_radtags` writes each individual sample to its own fastq file and produces a log file, all of which are located in `results/demultiplexed_fastqs/`. You can inspect this file by navigating to that directory and typing

```bash
less process_radtags.RADseq_tutorial_data.log
```

If you inspect individual fastq files, you should see that the barcode sequences have been removed, and each sequence begins with the SbfI recognition site. 

## Assess sequence quality with FastQC

`FastQC` is a program commonly used to assess the quality of Illumina sequence data. It reports on the base qualities, levels of adapter contamination, and various other features. 

Since we have 33 samples that can each be processed independently, we can submit this script as an array job. In an array job, we write a single script and the job scheduler executes it 33 times, substituting file names for each execution. When resources are available, this can significantly reduce time to completion. In this case it takes what would be a 10 minute job down to a minute or so. To understand how to specify array jobs, see our guidance [here](https://github.com/CBC-UCONN/Job-Arrays-on-Xanadu). 

The script is simple. The main call looks like this:

```bash
fastqc -t 2 -o $OUTDIR $INFILE
```

We can execute the script from the `scripts/lacewings` directory like this:

```bash
sbatch a2_fastqc.sh
```

When it is finished, we can use MultiQC to collate the FastQC output for each sample to facilitate checking the data:

```bash
multiqc ../../results/fastqc_dm
```

We execute the script as above:

```bash
sbatch a3_multiqc.sh
```

To inspect the results, we need to transfer the files to a local computer and view them in a web browser. We can use `scp` like this:

```bash
scp -r username@transfer.cam.uchc.edu:/full/path/to/RADseq_tutorial/results/fastqc_dm/multiqc* .
```

Or you can use a GUI-based FTP program like [FileZilla](https://filezilla-project.org/) or [Cyberduck](https://cyberduck.io/). 

On inspecting the data, we can see there are two samples with major problems. First, 107_downesi has substantial adapter contamination and a highly divergent GC content (these are all very closely related species, they should have very similar GC content!). We can keep it in the analysis, but if we look, we'll see later that it also shares very few variant loci with the other samples. The second problematic sample is 050_mediterranea. It has essentially dropped out of the analysis, having only 27 thousand sequence reads. 

We can carry these samples through most of the following steps, but they contribute much to the final analysis. 

