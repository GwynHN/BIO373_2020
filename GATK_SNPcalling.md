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

    $ module load Variants/GATK/3.8.1.0
    $ mkdir 03_callSNPs
    $ acc=516950
    $ java -jar $GATK_jar -T HaplotypeCaller -R 00_input/MedtrChr2.fa -I 02_dedup/${acc}.dedup.bam -ERC GVCF -nct 2 -o 03_callSNPs/${acc}.g.vcf.gz
    
This takes ~10 minutes for each sample. Maybe add this to your bash script and grab a coffee??

**Run again for sample 660389**


## Joint genotyping

This steps takes the input of gVCFs from each sample and combines them. This step recalculates some statistics and gives a more confident call of the variant at a particular site.

    $ java -jar $GATK_jar -T GenotypeGVCFs -R 00_input/MedtrChr2.fa -V 03_callSNPs/516950.g.vcf.gz -V 03_callSNPs/660389.g.vcf.gz -o 03_callSNPs/04_raw_variants.vcf.gz

# VCF format

All lines beginning with `##` contain information about how and when the VCF was generated and information about the flags included in the file. The single line with `#` tells you what each column of the following lines contains. Every line *without* a `#` is a variant.

The Genotype field (column 9) is important. Other flags may appear, but these are the minimum that should be included.\:

    
    GT: genotype 0/0 = homozygous ref; 0/1 = heterozygous; 1/1 = homozygous alt
    AD: allele depth; number of reads containing ref, alt base
    DP: total depth; total number of reads covering site
    GQ: genotype quality; difference between lowest and second lowest PL
    PL: genotype likelihood; what’s the probability it is NOT the correct genotype.
	The lowest is always 0.
    


:computer: You can view the contents of a VCF using zless. It is a regular text file that has been compressed. To search in the file using `zless`, type `/` when the file is loaded on the screen and type your search. i.e. `/CHROM` will take you to the line which shows the column names and the variants called in each sample. 



# Exercises Part 1

1. Here, we'll look in the VCF (04_raw_variants.vcf.gz) and take note of the information contained in the file (which is an overwhelming amount!). I like to get to the variants by searching for CHROM (`/CHROM`). You can look at any SNP, but I suggest searching for 7317, then 1018580. Those sites correspond to where we looked at in the BAM file in the mapping exercises. Take note of the variant quality (QD in INFO field). For an individual, take note of the genotype quality (GQ) and depth (AD and DP) as well. If you'd like, view the BAM file again using `samtools tview` and observe how the results we get from GATK compare to what you can see at those positions in a BAM file. Are the genotypes what you would expect just by looking at the BAM file?  

2. How many variants were discovered in this sample set?

3. Count the number of each genotype (0/0, 0/1, etc) for each sample. 


## Filtering out low quality data

This last step flags the variants that are low confidence. It will put the filterName in the INFO column and the G_filterName in the sample specific column if the site failed the filter. If at least one sample fails the G_filter, `FT` shows up in the FORMAT column. If all samples pass the G_filter, nothing shows up in the FORMAT column. This is something GATK might claim as a feature, not a bug. 

    $ java -jar $GATK_jar -T VariantFiltration -R 00_input/MedtrChr2.fa -o 03_callSNPs/05_variants_filtered.vcf.gz -V 03_callSNPs/04_raw_variants.vcf.gz --filterExpression "! vc.hasAttribute('QD') || QD < 2.0" --filterName "QD" --filterExpression "vc.isSNP() && (MQ < 30.0 || (vc.hasAttribute('MQRankSum') && MQRankSum < -15.0))" --filterName "MQ" -G_filter "GQ < 20 || DP == 0" -G_filterName "GQ" --setFilteredGtToNocall

# Exercises Part 2

1. Here, we'll look in the filtered VCF (05_variants_filtered.vcf.gz). This time, see if you notice what changed after the filtration step. For example, the FILTER field should now have a value (not just '.'). If you search for `FT`, you'll see in either samples the GT should have become a NoCall (./.).

2. Count the number of variants that failed each of the filters we applied (those applied to stats in the INFO field).

3. Count the number of genotypes that passed the genotype filter (GQ) for each sample. This one's a bit tricky due to the `FT` behavior noted above!

