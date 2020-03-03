# A few preliminaries. #

Before we get started with the tutorial, here are a few bits of background you should be aware of. 

- [Bash scripting](#bash-scripting)
  - [Pipes](#pipes)
  - [Shell variables](#shell-variables)
  - [`sed`](#sed)
- [The Xanadu cluster](#the-xanadu-cluster-and-slurm)
  - [What is Xanadu?](#what-is-xanadu)
  - [What is SLURM?](#what-is-slurm)
  - [The SLURM header](#the-slurm-header)
  - [Array jobs](#array-jobs)
- [Putting it all together](#putting-it-all-together)

## Bash scripting

This tutorial assumes you have some basic familiarity with `bash`. If you are not yet familiar with bash, see [here]() to learn and practice some of the basics. It also uses a few techniques that could be described as "intermediate bash scripting". Here are some of the more intermediate techniques we'll use in the tutorial:

### Pipes

bash scripts often chain together multiple commands with pipes: `|`. Pipes can be used with programs that can write output to the _standard output_ stream and read input from the _standard input_ stream. Try the following:

```bash
# echo writes to the standard output stream
echo Hello World!
# now redirect that output using a pipe to 'tr', which can read from the standard input and edit text
echo Hello World! | tr " " "\n"
# this turns the space into a line break. 
```
### Shell variables

We often use variables in shell scripting to make code easier to read. For example, in a command that looks like this:
```bash
	myprogram -I /very/extremely/long/path/to/GARBLEDACCESSIONX1848XKKD.fasta -O /another/very/extremely/long/path/to/RESULTSDIRECTORY/
```

We might clean it up with variables by:

```bash
INFILE=/very/extremely/long/path/to/GARBLEDACCESSIONX1848XKKD.fasta
OUTDIR=/another/very/extremely/long/path/to/RESULTSDIRECTORY/

myprogram -I $INFILE -O $OUTDIR
```
The shell variable is a text string and is interpreted in the context in which it is invoked. 

Note that when _setting_ variables you write `VARIABLE=something`, but when _invoking_ variables you write `$VARIABLE`. You need that `$` at the beginning. You may also invoke variables by ${VARIABLE}. 

Aside from making cleaner code, variables also provide the major advantage that we can easily substitute variables once in a script when necessary, instead of editing every occurence when a certain file, path or parameter needs to be changed. 

### `sed`
  
`sed`, or the "**s**tream **ed**itor" is a very flexible utility, but it is most often used to edit text read from the standard input stream. In this tutorial we use it to construct or edit input or output file names. 

```bash
# say we had an input file:
INFILE=007_downesi.fastq.gz
# and we wanted to create an output file name. 
# we can do:
OUTFILE=$(echo $INFILE | sed 's/.fastq.gz/.bam')
# now the variable $OUTFILE will contain "007_downesi.bam"
# and we could run our program
myprogram -I $INFILE -O $OUTFILE
```

Above, we define `$OUTFILE` by writing `$INFILE` to the standard output stream using `echo`, redirecting it to the standard input stream using `|`, and editing it using `sed`. We wrap the whole thing in `$()` to capture it and assign it to `$OUTFILE`. The find-replace operator for `sed` is `s/pattern/replacement/`. You can use regular expressions to do flexible pattern matching, but regexes are beyond the scope of this tutorial. 

## The Xanadu cluster and SLURM.

### What is Xanadu?

Xanadu is UConn's bioinformatics-oriented computer cluster. See our documentation on it [here](https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/xanadu/). 

### What is SLURM?

[SLURM](https://slurm.schedmd.com/documentation.html) is the job management software that coordinates requests for computational resources on Xanadu. Essentially, you write a script to execute some programs, submit it to SLURM, and SLURM sends it out to compute nodes in the cluster to be run. **ALL** work done on the cluster needs to be routed through SLURM. Again, see our documentation linked above. 

### The SLURM header. 

When you submit a script to SLURM, you need to add a header specifying the resources the work will need. We have detailed documentation [here](https://github.com/CBC-UCONN/CBC_Docs/wiki/Requesting-resource-allocations-in-SLURM), but briefly, a SLURM header looks like this:

```bash
#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=username@uconn.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-32]%20
```

The lines starting with `#SBATCH` tell `SLURM` what to do with the job and what resources it needs. 

### The `module` system. 

The system administrators of Xanadu install and maintain many pieces of software you might need, including `stacks`, and you can request new ones by [contacting us](https://bioinformatics.uconn.edu/). To manage all this software (including old versions), we use the `Environment Modules` package. To see a list of software available you can type `module avail`. This produces a long list. To use a given piece of software, you use `module load software/version`. To load `bwa`, for example you would type `module load bwa/0.7.17`. 

### Array jobs. 

A couple of scripts in this tutorial utilize _array jobs_ in SLURM. A job array is a series of very similar or identical jobs (called tasks). For this reason, they can be executed with a single script that may be modified slightly in each instance. A major advantage of array jobs is the ability to parallelize work without engaging in highly error-prone copy-paste-edit coding. 


## Putting it all together

To put all this together, we can use alignment of a collection of fastq files to a reference genome as an example. In the lacewing dataset we have 33 samples. To align them to the reference genome, we use [this script](/scripts/lacewings/c1_bwaalign.sh). It begins like this:

```bash
#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-32]%20

# print host name and date
hostname
date

# load software
module load bwa/0.7.17
module load samtools/1.9

# input, output files and directories
INDIR=../../results/demultiplexed_fastqs/

OUTDIR=../../results/aligned
mkdir -p $OUTDIR

# indexed reference genome
REFERENCE=/UCHC/PublicShare/CBC_Tutorials/RADseq_tutorial_data/lacewing_genome/redundans_metaquast_filtered.nomt.masked.fasta

# make a bash array of fastq files
FASTQS=($(ls -1 $INDIR/*.fq.gz))

# pull out a single fastq file
INFILE=$(echo ${FASTQS[$SLURM_ARRAY_TASK_ID]} | sed 's/.*\///')
# use sed to create an output file name
OUTFILE=$(echo $INFILE | sed 's/fq.gz/bam/')
```

First comes the SLURM header. The key to specifying this as an array job is the header line `#SBATCH --array=[0-32]%20`. This line means that the script should be run 33 times, with up to 20 tasks running at a time. Each task is given a number, from 0-32, which is assigned to the variable `$SLURM_ARRAY_TASK_ID`. We use each of the tools and concepts described above to load relevant software and define input and output files. There is one major wrinkle:

We use a _bash array_(`FASTQS`), which is essentially a list, to hold all the fastq file names. The array is specified like this:

`FASTQS=($(ls -1 $INDIR/*.fq.gz))`

In each instance of the script, we use `$SLURM_ARRAY_TASK_ID` to pull one element from the array and assign it to the variable `$INFILE`:

`INFILE=$(echo ${FASTQS[$SLURM_ARRAY_TASK_ID]} | sed 's/.*\///')`

The `sed` call uses a regular expression to remove the path from the file name. 

The actual call to the aligner is left out of this code snippet and will be described in the [relevant section](/01c_Lacewings_stacksrefmap.md). 

We also have written a [more detailed tutorial on using array jobs](https://github.com/CBC-UCONN/CBC_Docs/wiki/Job-arrays-on-Xanadu). 

