# Setting up for Variant calling with GATK

## Environment managememt with module system

Here, we'll again use the module system to load the software we want to use. You need to do this every time you log into the server.

    $ source /usr/local/ngseq/etc/lmod_profile
    
To check which modules (software) are available, type the following. Use arrow keys to navigate and type `q` to exit.

    $ module avail

## It's all relative

We will set up the directories in a specific way during the exercise and I have the commands set up to run based on this directory structure. I use relative paths a lot and most things will be run **assuming you are in the VarCall/ directory in YOUR directory on the server.** 

The most common reason things don't work are:

- Trying to run things from the wrong directory
- Unpaired quotation marks or brackets
- Misspellings 

Make sure you are in the directory you think you are in. Reminder: `pwd` command will tell you which directory you currently are in. 

Start from your directory on the server for the course and make a new folder for this exercise:

    $ ssh username@172.23.30.xx
    $ cd /scratch/bio373_2020/USERNAME
    $ mkdir VarCall

## Set up input by creating symlinks to input files

We use symlinks to help save diskspace on the server. 

    $ cd VarCall
    $ mkdir 00_input
    $ cd 00_input
    $ dataDir="/scratch/bio373_2020/data/SNPcalling/inputs"
    $ ln -s ${dataDir}/MedtrChr2.fa
    $ ln -s ${dataDir}/516950_chr2_R1.fastq.gz
    $ ln -s ${dataDir}/516950_chr2_R2.fastq.gz
    $ ln -s ${dataDir}/660389_chr2_R1.fastq.gz
    $ ln -s ${dataDir}/660389_chr2_R2.fastq.gz

### File types

**FASTA**

This is a file with nucleotide sequence(s) in it. For this exercise, it contains the reference sequence from part of chromosome 2 in *M. truncatula* to which we will map the reads. 

The line beginning with `>` is the information about the sequence that follows on the next line.

    $ less MedtrChr2.fa

More info: https://en.wikipedia.org/wiki/FASTA_format

**FASTQ** 

This is the output file given to you after sequencing and contain the reads.

    $ zless 516950_chr2_R1.fastq.gz

More info: https://en.wikipedia.org/wiki/FASTQ_format

For the quality score encoding: https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm

## Exercises

1. How many sequences are there in the reference sequence file (MedtrChr2.fa)? 

2. How many reads does each fastq file have (\*_R1.fastq.gz)? Does each sample have the same number of R1 and R2 reads? (Caution: Q scores can be + or @)

3. How many bases are in the reference sequence? How many missing bases (N)? Don’t forget ‘\n’ is considered a character!


## Answers 

On the server: `/scratch/bio373_2020/data/SNPcalling/answers/fastAQ.txt`
