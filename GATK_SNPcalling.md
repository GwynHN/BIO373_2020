# Using GATK to call variants
If you logged out of the server, you'll need to load the samtools and picard modules again.

    $ source /usr/local/ngseq/etc/lmod_profile
    $ module load Tools/samtools/1.10
    $ module load Tools/Picard/2.18.0


## Prep the reference sequence

GATK needs to use various index and dictionaries of the reference sequence in order to run smoothly. You need to use samtools and picard to generate these.

**Make sure you're in your VarCall/ directory**

    $ samtools faidx 00_input/MedtrChr2.fa
    $ java -jar $Picard_jar CreateSequenceDictionary R=00_input/MedtrChr2.fa O=00_input/MedtrChr2.dict

## Find some variants!

Now we use several tools in GATK to discover variants, group the samples, and finally filter the dataset to be used in downstream analyses.

First, you need to call variants in each samples. This creates a gVCF file. The step looks at the base that is called in the read that mapped at each site and determines whether there is a variant or not. 

   ## $ module load Variants/GATK/3.8.1.0
    $ module load Variants/GATK/4.1.8.0
    $ mkdir 03_callSNPs
    $ acc=516950
    $ java -jar $GATK_jar -T HaplotypeCaller -R 00_input/MedtrChr2.fa -I 02_dedup/${acc}.dedup.bam -ERC GVCF -nct 2 -o 03_callSNPs/${acc}.g.vcf.gz
    
    GATK4
% gatk HaplotypeCaller -R 00_input/MedtrChr2.fa -I 02_dedup/${acc}.dedup.bam -ERC GVCF -O 03_callSNPs/${acc}.g.vcf.gz
[September 21, 2020 6:08:25 PM CEST] org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done. Elapsed time: *12.06 minutes*. 
seems to use more than 1 thread...GATK4 has a ridiculous parallelization argument and takes longer!!!!!. what did they do?!?

**Run again for sample 660389**

Note: you can view the contents of a gVCF or VCF using zless. It is a regular text file that has been compressed. This is different than trying to view the BAM files before.

## Joint genotyping

This steps takes the input of gVCFs from each sample and combines them. This step recalculates some statistics and gives a more confident call of the variant at a particular site.

    $ java -jar $GATK_jar -T GenotypeGVCFs -R 00_input/MedtrChr2.fa -V 03_callSNPs/516950.g.vcf.gz -V 03_callSNPs/660389.g.vcf.gz -o 03_callSNPs/04_raw_variants.vcf.gz

# Exercises Part 1

1. Choose a SNP from the VCF (not gVCF) file (maybe one in the regions we looked at in the bam files earlier…). Take note of the variant quality (QD in INFO field). For an individual, take note of the genotype quality (GQ) and depth (AD and DP) as well. View the bam file using samtools tview and observe how the results we get from GATK compare to what you can see at those positions in a bam file. Are the genotypes what you would expect just by looking at the bam file? (7317 and 1018580 are searchable in the VCF)

2. How many variants were discovered in this sample set?

3. Count the number of each genotype (0/0, 0/1, etc) for each sample. 


## Filtering out low quality data

This last step flags the variants that are low confidence. It will put the filterName in the INFO column and the G_filterName in the sample specific column if the site failed the filter. If at least one sample fails the G_filter, `FT` shows up in the FORMAT column. If all samples pass the G_filter, nothing shows up in the FORMAT column. This is something GATK might claim as a feature, not a bug. 

    $ java -jar $GATK_jar -T VariantFiltration -R 00_input/MedtrChr2.fa -o 03_callSNPs/05_variants_filtered.vcf.gz -V 03_callSNPs/04_raw_variants.vcf.gz --filterExpression "! vc.hasAttribute('QD') || QD < 2.0" --filterName "QD" --filterExpression "vc.isSNP() && (MQ < 30.0 || (vc.hasAttribute('MQRankSum') && MQRankSum < -15.0))" --filterName "MQ" -G_filter "GQ < 20 || DP == 0" -G_filterName "GQ" --setFilteredGtToNocall

# Exercises Part 2

1. Count the number of variants that failed each of the filters we applied (those applied to stats in the INFO field).

2. Count the number of unique genotypes (0/0, 0/1, etc) that passed the genotype filter (GQ) for each sample.

