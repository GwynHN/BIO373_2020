# Mapping reads to reference genome: DNA edition

Here we map reads to a reference genome. These correspond to DNA, whereas Masa did RNA mapping/alignment with you last week. They are similar in the types of steps to execute, but the underlying concepts are a bit different.

Before you map, you usually want to check the quality of the reads using something like `fastqc` and then you definitely want to trim the reads. No use in trying to map poor quality data. These steps have already been done, so we start with the mapping  

We will use the module system to set up our environment:

    $ source /usr/local/ngseq/etc/lmod_profile

## Map reads with BWA

**Make sure you're back in your `VarCall/` directory**

First, set up the environment with the software we will use. BWA and samtools. In order to run the analysis, we need to index the reference sequence

    $ module load Aligner/BWA/0.7.17
    $ module load Tools/samtools/1.10
    $ bwa index 00_input/MedtrChr2.fa

We have two samples for this exercise, called `516950` and `660389`. Now we can run the mapping tool BWA and pipe the output through samtools, which will sort the output so we can use it in subsequent steps.

    $ mkdir 01_aligned
    $ acc=516950
    $ bwa mem -M -t 2 -R "@RG\tID:CAV90ANXX.6\tPL:Illumina\tLB:${acc}\tSM:${acc}" 00_input/MedtrChr2.fa 00_input/${acc}_chr2_R{1,2}.fastq.gz | samtools sort -m 16G -T /scratch/gwynhn -o 01_aligned/${acc}.sorted.bam
    $ samtools index 01_aligned/${acc}.sorted.bam

**Now change the value in ${acc} variable above for the other sample and re-run the whole bwa|samtools command**

:computer: Pro tip: if you press the up or down arrow keys, you can scroll through the commands you have already run during your shell session.

You should have 4 files in 01_aligned now: 2 `.bam` files and 2 `.bam.bai` files.

### BAM files and relative path practice

After mapping, you have a BAM file of where in the reference sequence the read mapped to. Because BAM files are compressed in a special way (not gzipped!) you can only open the BAM with certain tools, here `samtools`.

Below, I've given the basic commands but also a bit of thinking about where you are in the directory structure and how you can modify the input of the command so the computer knows where to look. 

Look at the header of the BAM file:

    $ samtools view -H 516950.sorted.bam 

(Are you in the correct directory? Can you modify the path to the file if not?)

Now look at the actual mapped reads:

    $ samtools view 516950.sorted.bam | less

Masa used IGV to view the alignments. You can also do that (remember you have to `scp` the BAM file onto the local computer!) Here's the command line alternative:

    $ samtools tview 516950.sorted.bam MedtrChr2.fa

(**Think of where the reference sequence is and where the bam file is RELATIVE to one another**)

While in the samtools tview, use the `?` to open the menu and `q` to exit. Once tview is open, you can type `g` and then within the prompt, type a location you'd like to go to. For example chr2:7271 then press `Enter`. `Esc` will get you out of the location search prompt. 

## Mark duplicates

    $ module load Tools/Picard/2.18.0
    $ mkdir 02_dedup
    $ java -jar $Picard_jar MarkDuplicates I=01_aligned/${acc}.sorted.bam O=02_dedup/${acc}.dedup.bam M=02_dedup/${acc}.metrics
    $ samtools index 02_dedup/${acc}.dedup.bam
    
**Run for the other sample!**

:computer: Forgot which sample number was stored in the acc variable? The echo command will show you!

    $ echo ${acc}

# Exercises

These exercises are really just designed to try to get you to understand what mapping does with the reads and what a potential variant might look like in a BAM file. There are fairly detailed.

Bitwise flag meaning: https://broadinstitute.github.io/picard/explain-flags.html

1. During mapping, each read is tagged with a bitwise flag that contains information about the read. Find bitwise flags for a few reads in any bam file and decode them using the link above. On the website, you can check only one box at a time to see what an individual property’s value is.

2. Using samtools tview, go to chr2:7271 in 516950.dedup.bam by pressing ‘g’ then typing chr2:7271[Enter]. Compare that with chr2:1018541. Why do you think they are different in terms of coverage and mapping quality? Find a few other areas that look interesting to you and take note of their position if you'd like to go back after the SNP calling step and see if a SNP was indeed called.

3. Do these regions look the same in sample 660389?

4. How many unique bitwise flags are there in 516950.sorted.bam file? The dedupped bam file? How many reads were marked as duplicates? (Hint: the flags are the sum of the value of each individual property assigned to a read; duplicate = 1024)

5. Write a bash script to run the alignment and dedup steps on both genotypes. To encourage organization and reproducibility, make a directory to keep your script(s) in and run from there :)

