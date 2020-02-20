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

1. Run `ustacks`. This step analyzes **each sample separately**. Reads are first sorted into groups with identical sequences. These groups fall into two classes based on a user-defined frequency threshold: "stacks", which it is hoped represent true allelic sequences, and "secondary reads", which it is hoped are mostly alleles with sequencing errors. If a targeted locus contains a heterozygous site, then at this stage it should be represented by two stacks. The allelic stacks are then assembled into putative loci if their sequence divergence falls below a user-specified threshold. The secondary reads are then mapped back to the loci, again contingent on a user-specified divergence threshold. 

2. Run `cstacks`. With a set of putative genetic loci identified **within** individuals we then run `cstacks` to merge loci across individuals and create a catalog of haplotypes for each locus. 

3. Run `sstacks`, which matches loci from individuals to the catalog created by `cstacks`. 

4. Run `tsv2bam`, which essentially sorts the data locus-wise. Even with paired-end sequencing data, all steps up to this point have used only the first read. In this step, the second reads are incorporated into bam files. At this stage we have clusters of sequences across individuals that correspond to putative genomic loci. 

5. Run `gstacks`. This module does three things: first, it assembles reads into contigs, which then form a kind of reference genome; second, it aligns the reads to those contigs; third, it calls variants and genotypes individuals against the reference contigs. 

6. Run `populations`, which calculates some basic statistics and filters and formats the data for downstream analysis.

For further details about the implementation of `Stacks`, see Rochette et al. (2017) an Catchen et al. (2011). 

It is worth noting that there are many parameters users may set to influence the assembly of loci and calling of variants and it is not always clear _a priori_ what they should be set to. Paris et al. (2017) recommend searching across parameter space to find the set of parameters that maximize the number of polymorphic loci found in >=80% of individuals. 

## Step 1: ustacks

`ustacks` clusters sequences and identifies loci within each sample. This means we can parallelize it by executing it separately for each sample. This is an advantage of running the separate modules on their own. 

Here we'll run `ustacks` using [this script](/scripts/lacewings/a4_ustacks.sh). Before explaining how to parallelize it, we'll first look at the `ustacks` call and its options. 

```bash
ustacks \
-f $INFILE \
-o $OUTDIR \
-i $ID \
--name $SAM \
-t gzfastq \
-p 6 \
-M 8 \
-m 3 \
--max-gaps 10 \
--high-cov-thres 10 
```

In this code, some key files are represented as shell variables (e.g. `$INFILE`) that are defined higher in the script. 

The parameters dealing with file management and program execution are as follows: `-f` gives the input fastq file, `-o` gives the output directory, `-i` is an integer, which should be unique to each sample, `--name` gives the sample ID, `-t` indicates input fastq files are gzipped amd `-p` specifies the number of cpu threads `ustacks` should use. 

Several parameters may be modified to change how the algorithm assembles loci within each sample: `-M` is the number of nucleotide differences above which allelic stacks will not be merged into loci. `-m` is the minimum number of reads required to create a stack, and `-N`, which we do not specify, and so leave at the default, is the maximum number of differences allowed to assign a "secondary read" to a stack. 

Because this dataset has very high coverage and multiple species, some with high genetic diversity, we have raised `-m` and `-M` above their default values. 

## Step 2: cstacks


## Step 3: sstacks


## Step 4: tsv2bam


## Step 5: gstacks


## Step 6: populations


## References

Catchen, Julian M., Angel Amores, Paul Hohenlohe, William Cresko, and John H. Postlethwait. 2011. “Stacks: Building and Genotyping Loci de Novo from Short-Read Sequences.” G3  1 (3): 171–82.

Paris, Josephine R., Jamie R. Stevens, and Julian M. Catchen. 2017. “Lost in Parameter Space: A Road Map for Stacks.” Edited by Susan Johnston. Methods in Ecology and Evolution / British Ecological Society 8 (10): 1360–73.

Rochette, Nicolas C., Angel G. Rivera-Colón, and Julian M. Catchen. 2019. “Stacks 2: Analytical Methods for Paired-End Sequencing Improve RADseq-Based Population Genomics.” Molecular Ecology 28 (21): 4737–54.

