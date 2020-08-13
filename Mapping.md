# Mapping reads to reference genome: DNA edition

Here we map reads to a reference genome. These correspond to DNA, whereas Masa did RNA mapping/alignment with you last week. They are similar in the types of steps to execute, but the underlying concepts are a bit different.

Before you map, you usually want to check the quality of the reads using something like `fastqc` and then you definitely want to trim the reads. No use in trying to map poor quality data. These steps have already been done, so we start with the mapping  

In case:

    $ source /usr/local/ngseq/etc/lmod_profile

## Map reads with BWA

**Make sure you're back in your VarCall/ directory**

First, set up the environment with the software we will use. BWA and samtools. In order to run the analysis, we need to index the reference sequence

    $ module load Aligner/BWA/0.7.15
    $ module load Tools/samtools/1.9
    $ bwa index 00_input/MedtrChr2.fa

Now we can run the mapping tool BWA and pipe the output through samtools, which will sort the output so we can use it in subsequent steps. You have two samples to process: 516950 and 660389.

    $ mkdir 01_aligned
    $ acc=516950
    $ bwa mem -M -t 2 -R "@RG\tID:CAV90ANXX.6\tPL:Illumina\tLB:${acc}\tSM:${acc}" 00_input/MedtrChr2.fa 00_input/${acc}_chr2_R{1,2}.fastq.gz | samtools sort -m 16G -T /scratch/gwynhn -o 01_aligned/${acc}.sorted.bam
    $ samtools index 01_aligned/${acc}.sorted.bam

**Now change the value in ${acc} variable above and re-run the whole bwa|samtools command**

Don't forget you can use `echo` to see what you have stored in the `acc` variable. 

You should have 4 files in 01_aligned now: 2 .bam files and 2 .bam.bai files.

### BAM files and relative path practice

After mapping, you have a BAM file of where in the reference sequence the read mapped to. Because BAM files are compressed in a special way (not gzipped!) you can only open the BAM with certain tools, here `samtools`.

Below, I've given the basic commands but also a bit of thinking about where you are in the directory structure and how you can modify the input of the command so the computer knows where to look. 

Look at the header of the BAM file:

    $ samtools view -H 516950.sorted.bam 

(Are you in the correct directory? Can you modify the path to the file if not?)

Now look at the read in the BAM file and how the mapping is :
    
    $ samtools view 516950.sorted.bam | less
    # type `q` to get out of less view mode
    
Masa used IGV to view the alignments. You can also do that (remember you have to scp the BAM file onto the local computer!) Here's the command line alternative:

    $ samtools tview 01_aligned/516950.sorted.bam 00_input/MedtrChr2.fa
    
The top line is the position along the reference sequence, the 2nd line is the reference sequences (MedtrChr2), the 3rd line is the consensus sequence based on the reads, and the rest of the lines are the reads and where they mapped.

While in the samtools tview, use the `?` to open the menu and `q` to exit. Once tview is open, you can type `g` and then within the prompt, type a location you'd like to go to. For example chr2:7271 then press `Enter`. `Esc` will get you out of the location search prompt. 

## Mark duplicates

    $ module load Tools/Picard/2.18.0
    $ mkdir 02_dedup
    $ java -jar $Picard_jar MarkDuplicates I=01_aligned/516950.sorted.bam O=02_dedup/516950.dedup.bam M=02_dedup/516950.metrics
    $ samtools index 02_dedup/516950.dedup.bam
**Run for the other sample!**

# Exercises

These exercises are really just designed to try to get you to understand what mapping does with the reads and what a potential variant might look like in a BAM file. This is not part of the exam.

Bitwise flag meaning: https://broadinstitute.github.io/picard/explain-flags.html

1. Find bitwise flags for a few reads in any bam file and decode them using the link above. On the website, you can check only one box at a time to see what an individual property’s value is.

2. Using samtools tview, go to chr2:7271 in 516950.dedup.bam by pressing ‘g’ then typing `chr2:7271[Enter]`. Compare that with chr2:1018541. Why do you think they are different in terms of coverage and mapping quality? Find a few other areas that look interesting to you and take note of their position.

3. Do these regions look the same in sample 660389?

4. How many unique bitwise flags are there in 516950.sorted.bam file? The dedupped bam file? How many reads were marked as duplicates? (Hint: the flags are the sum of the value of each individual property assigned to a read; duplicate = 1024)

5. Write a script to run the alignment and dedup steps on both genotypes.

