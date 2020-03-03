# A few preliminaries. #

Before we get started with the tutorial, here are a few bits of background you should be aware of. 

## Bash scripting

This tutorial assumes you have some basic familiarity with `bash` and uses a few techniques that could be described as "intermediate bash scripting". If you are not yet familiar with bash, see [here]() learn and practice some of the basics. Here are some of the more intermediate techniques we'll use in the tutorial:

- Pipes

  bash scripts often chain together multiple commands with pipes: `|`. Pipes can be used with programs that can write output to the _standard output_ stream and read input from the _standard input_ stream. Try the following:

  ```bash
  # echo writes to the standard output stream
  echo Hello World!
  # now redirect that output using a pipe to 'tr', which can read from the standard input and edit text
  echo Hello World! | tr " " "\n"
  # this turns the space into a line break. 
  ```
- Shell variables

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

  Note that when _setting_ variables you write `VARIABLE=something`, but when _invoking_ variables you need that `$` at the beginning. You may also invoke variables by ${VARIABLE}. 

  Aside from making cleaner code, variables also provide the major advantage that we can easily substitute variables once in a script when necessary, instead of editing every occurence when a certain file, path or parameter needs to be changed. 

- `sed`
  
  `sed`, or "the **s**tream **ed**itor" is a very flexible utility, but it is most often used to edit text read from the standard input stream. In this tutorial we might use it to construct or edit input or output file names. 

  ```bash
  # say we had an input file:
  INFILE=007_downesi.fastq.gz
  # and we wanted to create an output file name. 
  # we can do:
  OUTFILE=$(echo $INFILE | sed 's/.fastq.gz/.bam')
  # now the variable outfile will contain "007_downesi.fastq.gz"
  ```
  Above, we write the input file name to the standard output stream using echo, redirect it to the standard input stream using a pipe, and edit it using `sed`. The find-replace operator for `sed` is `s/pattern/replacement/`. You can use regular expressions to do flexible pattern matching, but regexes are beyond the scope of this tutorial. 

xanadu and SLURM
- what is xanadu?
- what is SLURM?
- SLURM header
- "modules"

array jobs
- what is an array job and how does it work?

A couple of scripts in this tutorial utilize _array jobs_ in SLURM. 