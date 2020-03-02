# Running Stacks in _de novo_ mode. #

## Introduction

This section of the tutorial deals with running `Stacks` in _de novo_ mode on the lacewing dataset. 

Steps here will use the following software packages:

- [ Stacks ](http://catchenlab.life.illinois.edu/stacks/)
- [ Stacks ](http://catchenlab.life.illinois.edu/stacks/)

Each major step has an associated bash script tailored to the UConn CBC Xanadu cluster with appropriate headers for the [Slurm](https://slurm.schedmd.com/documentation.html) job scheduler. The code can easily be modified to run interactively, or using another job scheduler. 

## Contents
  
1.    [ Motivation ](#Motivation)
2.    [ Step 1: ustacks ](#Step-1-ustacks)
3.    [ Step 2: cstacks ](#Step-2-cstacks)
4.    [ Step 3: sstacks ](#Step-3-sstacks)
5.    [ Step 4: tsv2bam ](#Step-4-tsv2bam)
6.    [ Step 5: gstacks ](#Step-5-gstacks)
7.    [ Step 6: populations ](#Step-6-populations)

## Motivation

`Stacks` implements two methods for detecting variants in sets of sequenced samples: a _de novo_ approach, in which loci are assembled by sequence similarity, and a reference-mapping approach, in which loci are identified after alignment of all sequences to a reference genome. Here we will cover the _de novo_ approach, and demonstrate the reference mapping approach in the next section. 

`Stacks` implements the entire _de novo_ pipeline in a single perl script `denovo_map.pl`, but here we will run each constituent module separately to demonstrate it. 

Before we get started, it's worth explaining a little bit about how `Stacks` works in _de novo_ mode. There are six submodules to run. The first four collectively identify sequences within and across individuals that correspond to homologous loci, the fifth and sixth generate, filter and output genotypes. They are as follows:

1. `ustacks`: This step analyzes **each sample separately**. Reads are first sorted into groups with identical sequences. These groups fall into two classes based on a user-defined frequency threshold: "stacks", which it is hoped represent true allelic sequences, and "secondary reads", which it is hoped are mostly alleles with sequencing errors. If a targeted locus contains a heterozygous site, then at this stage it should be represented by two stacks. The allelic stacks are then assembled into putative loci if their sequence divergence falls below a user-specified threshold. The secondary reads are then mapped back to the loci, again contingent on a user-specified divergence threshold. 

2. `cstacks`: With a set of putative genetic loci identified **within** individuals we then run `cstacks` to merge loci across individuals and create a catalog of haplotypes for each locus. 

3. `sstacks`: This module matches loci from individuals to the catalog created by `cstacks`. 

4. `tsv2bam`: This module essentially sorts the data locus-wise instead of sample-wise. Even with paired-end sequencing data, all steps up to this point have used only the first read. In this step, the second reads are incorporated into bam files. At this stage we have clusters of sequences across individuals that correspond to putative genomic loci. 

5. `gstacks`: This module does three things: first, it assembles reads into contigs, which then form a kind of reference genome; second, it aligns the reads to those contigs; third, it calls variants and genotypes individuals against the reference contigs. 

6. `populations`: This step calculates some basic statistics and filters and formats the data for downstream analysis.

For further details about the implementation of `Stacks`, see Rochette et al. (2017) an Catchen et al. (2011). 

It is worth noting that there are many parameters users may set to influence the assembly of loci and calling of variants and it is not always clear _a priori_ what they should be set to. Paris et al. (2017) recommend searching across parameter space to find the set of parameters that maximize the number of polymorphic loci found in >=80% of individuals. 

## Step 1: ustacks

`ustacks` clusters sequences and identifies loci within each sample. This means we can parallelize it by executing it separately for each sample. This is an advantage of running the separate modules on their own. 

Here we'll run `ustacks` using [this script](/scripts/lacewings/b1_ustacks.sh). Before explaining how to parallelize it, we'll first look at the `ustacks` call and its options. 

```bash
ustacks \
-f $INFILE \
-o $OUTDIR \
-i $ID \
--name $SAM \
-t gzfastq \
-p 6 \
-M 6 \
-m 3 \
--max-gaps 5 \
--high-cov-thres 10 
```

In this code, some key files are represented as shell variables (e.g. `$INFILE`) that are defined higher in the script. 

__Parameters dealing with file management and program execution__: `-f` gives the input fastq file, `-o` gives the output directory, `-i` is an integer, which should be unique to each sample, `--name` gives the sample ID, `-t` indicates the format of the input files (here gzipped fastq) and `-p` specifies the number of cpu threads `ustacks` should use. `-p` should be the same as the number you set in the SLURM header: `#SBATCH --cpus-per-task=6`). 

__Algorithmic parameters impacting assembly__: `-M` is the number of nucleotide differences above which allelic stacks will not be merged into loci. `-m` is the minimum number of reads required to create a stack, and `-N`, which we do not specify, and so leave at the default, is the maximum number of differences allowed to assign a "secondary read" to a stack. 

You can think of the `-M` parameter as being related to expected heterozygosity. This number varies widely among species (Leffler et al. 2012). Humans are at about 0.001, while some marine invertebrates can range up to 0.05. A good starting place for this parameter might be 2 times the expected heterozygosity times the read length. Setting the value too high will cause paralogous loci to be collapsed, and setting it too low will cause homologous loci to be split. Because this dataset has very high coverage and multiple species, some with high genetic diversity, we have raised `-m` and `-M` above their default values. 

The parameter `--max-gaps` gives the maximum number of gaps that can exist between two stacks to merge them. We set this higher than the default here because these are 150bp reads. `--high-cov-thre` the number of standard deviations above the mean that a stack can be before it is discarded as a repetitive element. We set this higher for this dataset because preliminary analyses at 3 (the default) were discarding 80% of the data. 

`ustacks` will create three files for each sample: `samplename.alleles.tsv.gz`, `samplename.snps.tsv.gz`, and `samplename.tags.tsv.gz`. All will be located in the single specified output directory (see their contents [here](http://catchenlab.life.illinois.edu/stacks/manual/#files)). 

__Parallelizing ustacks__: In this script, we parallelize `ustacks` using a feature of the SLURM job scheduler called _job arrays_. For a detailed explanation of how to use job arrays, see [here](https://github.com/CBC-UCONN/CBC_Docs/wiki/Job-arrays-on-Xanadu). Briefly, in a job array we specify an extra line in the SLURM header:

`#SBATCH --array=[0-32]%20`

This indicates the script should be run 33 times, with 20 jobs allowed at once. In each instance, the shell variable `SLURM_ARRAY_TASK_ID` will be given a number between 0 and 32. We use these numbers to pull out a different sample and assign it to the `SAM` variable and to set the `ID` variable. This is advanced cluster usage, and you don't need to know how to do it, but it can be very helpful. If you **don't** wish to run this step in parallel, you can run `denovo_map.pl` giving it a list of all your samples and you won't have to deal with job arrays. 

We can run the script by navigating to the `scripts/lacewings/` directory and typing

```bash
sbatch b1_ustacks.sh
```

One of the important things to have a look at here are the summaries of how the data was assembled. In our script, these summaries are written to slurm `.err` files. There is one for each task in the job array, and they'll be located in the directory `scripts/lacewings`. The actual job numbers will change, but here we can have a look at `ustacks_724336_1.err`. This should correspond to the sample `010_downesi`. The file first lists the parameters that were run and then output that looks like this:

```
Loaded 2679865 reads; formed:
  18944 stacks representing 2421459 primary reads (90.4%)
  241086 secondary stacks representing 258406 secondary reads (9.6%)

Stack coverage: mean=127.82; stdev=1735.81; max=152812; n_reads=2421459(90.4%)
Removing repetitive stacks: cov > 17486 (mean+10*stdev)...
  Blacklisted 2740 stacks.
Coverage after repeat removal: mean=102.25; stdev=264.57; max=17475; n_reads=1656879(61.8%)

Assembling stacks (max. dist. M=8)...
  Assembled 16204 stacks into 7914; blacklisted 314 stacks.
Coverage after assembling stacks: mean=155.85; stdev=199.45; max=1601; n_reads=1182929(44.1%)

Merging secondary stacks (max. dist. N=10 from consensus)...
  Merged 180024 out of 258406 secondary reads (69.7%), 60738 merged with gapped alignments.
Coverage after merging secondary stacks: mean=179.57; stdev=240.74; max=8373; n_reads=1362953(50.9%)

Assembling stacks, allowing for gaps (min. match length 80.0%)...
  Assembled 7914 stacks into 6727 stacks.
Coverage after gapped assembly: mean=212.86; stdev=261.43; max=8373; n_reads=1362953(50.9%)

Merging secondary stacks, allowing for gaps (min. match length 80.0%)...
  Merged 63221 out of 78382 secondary reads (80.7%).
Coverage after merging gapped secondary stacks: mean=222.74; stdev=277.72; max=8373; n_reads=1426174(53.2%)

Final coverage: mean=222.74; stdev=277.72; max=8373; n_reads=1426174(53.2%)
```

There are six steps here. 

1. Creating the initial stacks (90.4% of reads) and secondary read stacks (9.6% of reads). 
2. Removing excessively high coverage stacks. This throws away 28.6% of the data! You can see that the average initial stack depth for **this** dataset is 102.25x. `-m` is three, but could reasonably be set **much** higher. 50 would probably be reasonable. This is extremely dataset dependent. 
3. Assembling the stacks into loci and again removing high coverage stacks. Now we've lost another 17.7% of the data!
4. Merging secondary stacks. 
5. Merging stacks with gaps. 
6. Merging secondary stacks with gaps. 

In the end, only 53% of the data are assembled into 6727 potentially usable loci. The rest of the samples are consistent with this. I have tried to modify the parameters to do better, but the variance in coverage among stacks is extreme. This is either due to a highly repetitive genome or issues with sample preparation. With single-end RAD we don't have a lot of ways to diagnose this. 

## Step 2: cstacks

`cstacks` now takes all the putative genetic loci identified within individuals, compares them against each other and merges them to create a catalog of loci across species. 

This step cannot be parallelized, as it needs to take in information across all samples. Here we'll run `cstacks` using [this script](/scripts/lacewings/b2_cstacks.sh). The call looks like this:

```bash
cstacks \
-P $INDIR \
-M ../../metadata/lacewing_popmap.txt \
-p 20 \
--max-gaps 10 \
-n 15
```

This step has many fewer options. Again, we specify the input directory as a shell variable (defined higher in the script) `-P $INDIR`. Following that, we have the population map, `-M`. This is a tab-separated file with two columns: first, the sample names and second, the populations they belong to. In our case, the first few lines look like this:

```
135_downesi	downesi
007_downesi	downesi
010_downesi	downesi
160_downesi	downesi
161_downesi	downesi
114_carnK	carnK
111_carnK	carnK
112_carnK	carnK
109_carnK	carnK
113_carnK	carnK
...
```

Then we specify the number of CPU threads to use with `-p`. Here we've specified `--max-gaps 10` to allow loci with many or larger gaps to be merged. We are trying to make the assembly parameters more lax to deal with a multi-species dataset. Finally, we set `-n 15`. This flag is the maximum allowable number of differences between loci from different samples. The `Stacks` documentation suggests this number be identical to `-M` in the previous step, but as we expect genetic divergence between species to be higher than allelic divergence within, we set this higher here, to the equivalent of 10% divergence over a 150bp read. 

We can run the script by navigating to the `scripts/lacewings/` directory and typing

```bash
sbatch b2_cstacks.sh
```

## Step 3: sstacks

`sstacks` now matches loci from each sample to the catalog created by `cstacks`. This step can be parallelized as in `ustacks`, but it is much faster, so unless you have many samples, it may not be worth it. We'll run `sstacks` using [this script](/scripts/lacewings/b3_sstacks.sh). It's a simple call:

```bash
sstacks -P $INDIR -M $POPMAP -p 20
```
It uses the same population map as above, here specified as a shell variable, and we are setting `-p 20` to allow it to use 20 CPU threads. 

We can run the script by navigating to the `scripts/lacewings/` directory and typing

```bash
sbatch b3_sstacks.sh
```

## Step 4: tsv2bam

`tsv2bam` now sorts the data for each sample by locus, so variants can be jointly called across samples more efficiently. We'll run `tsv2bam` using [this script](/scripts/lacewings/b4_tsv2bam.sh). The call looks like this:

```bash
tsv2bam -P $INDIR -M $POPMAP -t 20
```
The options are the same as above, except `-t` gives the number of CPU threads. 

We can run the script by navigating to the `scripts/lacewings/` directory and typing

```bash
sbatch b4_tsv2bam.sh
```

## Step 5: gstacks

`gstacks` now assembles loci, aligns reads against them, and calls variable sites and individual genotypes. We'll run `gstacks` using [this script](/scripts/lacewings/b5_gstacks.sh). 

```bash
gstacks -P $INDIR -M $POPMAP -t 20
```
Here we use the same options as above, but some can be altered to influence how variant calling is done. 

We can run the script by navigating to the `scripts/lacewings/` directory and typing

```bash
sbatch b5_gstacks.sh
```
## Step 6: populations

Finally, we can filter and reformat the data for use in downstream applications using `populations`. We'll run `populations` using [this script](/scripts/lacewings/b6_populations.sh). The call looks like this:

```bash
populations \
-P $INDIR \
-M $POPMAP \
-p 5 \
-r 2 \
--hwe \
--genepop \
--vcf \
--treemix \
-t 8
```

The input and population files are specified as above. `-p` sets the minimum number of populations a locus must appear in to be output. `-r` is the minimum number of samples per population a locus must appear in to be output. `--hwe` calculates the departure from Hardy-Weinberg equilibrium for each locus. The next three options output the genotypes in three different formats and `-t` indicates the number of CPU threads that should be used. 

We can run the script by navigating to the `scripts/lacewings/` directory and typing

```bash
sbatch b6_populations.sh
```
This should be pretty quick, and at the end we should have a number of output files that all start with "populations" in the `results/stacks/denovo` directory. 


## References

Catchen, Julian M., Angel Amores, Paul Hohenlohe, William Cresko, and John H. Postlethwait. 2011. “Stacks: Building and Genotyping Loci de Novo from Short-Read Sequences.” G3  1 (3): 171–82.

Leffler, Ellen M., et al. "Revisiting an old riddle: what determines genetic diversity levels within species?." PLoS biology 10.9 (2012).

Paris, Josephine R., Jamie R. Stevens, and Julian M. Catchen. 2017. “Lost in Parameter Space: A Road Map for Stacks.” Edited by Susan Johnston. Methods in Ecology and Evolution / British Ecological Society 8 (10): 1360–73.

Rochette, Nicolas C., Angel G. Rivera-Colón, and Julian M. Catchen. 2019. “Stacks 2: Analytical Methods for Paired-End Sequencing Improve RADseq-Based Population Genomics.” Molecular Ecology 28 (21): 4737–54.

