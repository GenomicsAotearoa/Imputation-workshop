# Aim
1. Go through the pipeline of phasing and imputing high-density genotype to sequence level
2. Understand the importance of quality control 
3. Know how to evaluate the imputation performance

# Tools we need

>All of these tools/applications are pre-installed in the training platform 

BCFtools: basic bioinformatics software, in this tutorial, we use it for creating subsets and quality control. the current version is 1.10.2-GCC-9.2.0 ([http://samtools.github.io/bcftools/bcftools.html](http://samtools.github.io/bcftools/bcftools.html))

VCFtools: basic bioinformatics software, in this tutorial, we use it for comparing two vcf files and evaluate the concordance rate. The current version is /0.1.15-GCC-9.2.0-Perl-5.30.1 ([https://vcftools.github.io/index.html](https://vcftools.github.io/index.html))

R: basic statistics software

Beagle: software for phasing and imputation. The current version is version 5.1-18May20. ([https://faculty.washington.edu/browning/beagle/beagle.html](https://faculty.washington.edu/browning/beagle/beagle.html))

Minimac3: software for imputation ([https://genome.sph.umich.edu/wiki/Minimac3](https://genome.sph.umich.edu/wiki/Minimac3))

# Input files
In this tutorial, we used 5MB segment (30-35MB) of chr13 from the 1000 Genome phase 3, which is publicly available and can be downloaded from https://mathgen.stats.ox.ac.uk/impute/impute_v2.html

# Procedures
## 1. Load the module that required on NeSI
To find modules which are already available on NeSI, use `module spider #module_name`

To load modules and start using, use `module load #module_name`

``` bash
module load BCFtools/1.10.2-GCC-9.2.0
module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1
module load R/4.0.1-gimkl-2020a
module load Beagle/5.1-18May20.d20 
#minimac3=/nesi/nobackup/nesi02659/SEP28/practical2/Minimac3
```

## 2. Copy the files into the home directory
Define two directories: workshop directory and home directory. In this workshop, the analysis will be conducted in the home directory

```
maindir=/nesi/nobackup/nesi02659/SEP28/practical2
cd $HOME
```
```
mkdir -p imputation_workshop
cd $HOME/imputation-workshop
```

The main input file (seqvcf) is extracted (from 30MB to 35MB) on chr13 from 1000 human genome data. I selected some of the reliable SNPs to generate a HD genotype dataset(hdvcf). 

```
seqvcf=$maindir/nonfilter_seq_5MB.vcf.gz
hdvcf=$maindir/hd_5MB.vcf.gz
```

Now what you need to do is to cp both genotype and sequence data to your own home directory. In addition, we have to software that we need to use which are not available on NeSI. Please also copy these two files in your own working directory.

```
cp $seqvcf $HOME/imputation-workshop
cp $hdvcf $HOME/imputation-workshop
cp $beagle5 $HOME/imputation-workshop
cp $minimac3 $HOME/imputation-workshop
```

The genotype and sequence files use "vcf.gz" format. We can not open it directly. To check how the genotype looks like, you need to use: zless -S. -S is to make the file well formated. 

```
zless -S $seqvcf
zless -S $hdvcf
```

## 3. Explore the input file

To have a basic idea of the genotype, BCFtools have a very convenient function: `stats`. By checking the original sequence file's information. the code you need is as below. -s is a common tag to show "samples". Samples to include or "-" to apply all variants. Via adding this, we also generate the statistics for each individual. The output will be named: original.stats

`bcftools stats -s "-" nonfilter_seq_5MB.vcf.gz > original.stats`

Now let's have a look at the output:

`less original.stats`

![](https://github.com/GenomicsAotearoa/Imputation-workshop/blob/master/Tutorial/img/Screen%20Shot%202020-09-16%20at%2014.37.07.png?raw=true)

So in the first part, you can see the basic statistics of the sequence file. We have 2504 samples in total, with 149,854 variants. 144,069 of them are SNPs, 5777 indels, 98 other, 694 multi-allelic sites with 467 multi-allelic SNPs.  

## 4. Quality control process

There are a lot of parameters you may take into consideration in your dataset, such as non-variant, singletons, multi-allelic positions, map quality, mendelian error, minor allele frequency, etc. The parameter settings vary depends highly on data, and also for the purpose of data. For example, if the data is for GWAS analysis, you probably don't want to throw away too many variants since a lot of the causals are rare alleles. If your data is used for prediction, then you probably don't want to include too much crap in the data, which is not beneficial for the next step model construction.

```
bcftools view -O z -o nosingleton.vcf.gz -i 'AC>1' nonfilter_seq_5MB.vcf.gz
tabix -f nosingleton.vcf.gz
bcftools stats -s "-" nosingleton.vcf.gz > step1.stats
less step1.stats
```

![](https://github.com/GenomicsAotearoa/Imputation-workshop/blob/master/Tutorial/img/Screen%20Shot%202020-09-16%20at%2014.38.50.png?raw=true)

As we can see here, after removing the non-variants and singletons, the number of variants decreased to 85,928, 80,184 SNPs, 5767 indels, 67 others, 694 multiallelic sites with 467 multi-allelic SNPs.  

```
bcftools view -O z -o nosingleton_2alleles.vcf.gz --max-alleles 2 nosingleton.vcf.gz
tabix -f nosingleton_2alleles.vcf.gz
bcftools stats -s "-" nosingleton_2alleles.vcf.gz > step2.stats
less step2.stats
```
![](https://github.com/GenomicsAotearoa/Imputation-workshop/blob/master/Tutorial/img/Screen%20Shot%202020-09-16%20at%2014.40.26.png?raw=true)

As we can see here, after setting the maximum allele into 2, the multi-allelic variants should be gone. The total number of variants decreased to 85,234, 79,627 SNPs, 5543 indels and 64 others. In this tutorial, I am only gonna consider non-variant, singleton and multi-allelic variants since the majority of the imputation software can not handle them anyway. Other filtering processes can also be done using bcftools, you can just pop in the website and go to filtering sessions. 

## 5. Generate the data for the imputation process

We are gonna just use this data for imputation. So there are two populations involved: reference population and study population. In reality you should have both datasets ready. In this tutorial, I decide to just use this one dataset. Treat the first 1000 samples as my reference, and the rest will be set as my study population. Using bcftools to extract the sample ID and basic awk function to generate two population ID files as follow:

```
bcftools query -l nonfilter_seq_5MB.vcf.gz > seq.ID
awk 'NR>0&&NR<=1000' seq.ID > ref.ID
awk 'NR>1000&&NR<=2000' seq.ID > study.ID
```

Now we create a new directory for our imputation:

```
mkdir -p imputation
cd imputation
```

bcftools is very handy of extracting the sample from the whole dataset via using -S. We will generate two files for our study samples: one from HD genotype as my study population. And another one from filtered sequence data for future validation. After extraction, we use function tabix to index both files.

```
bcftools view -O z -o study_hd.vcf.gz -S $HOME/study.ID $HOME/hd_5MB.vcf.gz
tabix -f study_hd.vcf.gz
bcftools view -O z -o study_filtered.vcf.gz -S $HOME/study.ID $HOME/nosingleton_2alleles.vcf.gz #for validation
tabix -f study_filtered.vcf.gz
```

Similarly, we will extract the samples for our reference samples. The same functions will be used here. We will generate two files for the reference samples, one is from the original unfiltered sequence, and another one from the filtered sequence. Again, we will use tabix to index these two files. 

```
bcftools view -O z -o ref_nonfiltered.vcf.gz -S $HOME/ref.ID $HOME/nonfilter_seq_5MB.vcf.gz
tabix -f ref_nonfiltered.vcf.gz
bcftools view -O z -o ref_filtered.vcf.gz -S $HOME/ref.ID $HOME/nosingleton_2alleles.vcf.gz
tabix -f ref_filtered.vcf.gz
```

## 6. Phasing using Beagle 5 

For imputation, no matter which software are you using, phasing is compulsory for the reference population. For some software like minimac3, both study and reference population have to be phased.  If you simply pop the unphased reference population into your imputation software, the software will immediately give you an error message. 

The data I downloaded already finished phasing that you can see in the dataset, all the genotypes are phased (eg: "0|1", not "0/1"). Besides phasing needs a lot of computation resources and a certain amount of time. Here since the data is already phased, it won't take too long.  

```
beagle gt=ref_filtered.vcf.gz chrom=13 out=ref_filtered_phased
tabix -f ref_filtered_phased.vcf.gz
```

```
beagle gt=ref_nonfiltered.vcf.gz chrom=13 out=ref_nonfiltered_phased
tabix -f ref_nonfiltered_phased.vcf.gz
```

## 7. Imputation using Beagle 5

In this tutorial, I will show you the imputation using two software: Beagle 5 and minimac3. Both software are very stable, reliable, easy-to-use, free, and pretty popular. Beagle 5 is computationally demanding but can give you accurate results very fast. Minimac is computationally efficient, but a bit slower. In addition, to prove that quality control is an important procedure, both filtered reference and unfiltered sequence reference will be used.

Beagle has been evolved from version 3.0 to the current 5.1 version. It becomes much faster and simpler. And be able to handle large datasets. In th e meantime, the parameters for running the software have been reduced a lot. There are several important parameters that can influence imputation performance such as effective population size (Ne), window size, etc. Check the following paper: Improving Imputation Quality in BEAGLE for Crop and Livestock Data [https://www.g3journal.org/content/10/1/177](https://www.g3journal.org/content/10/1/177)

```
beagle gt=study_hd.vcf.gz ref=ref_filtered_phased.vcf.gz chrom=13 impute=true gp=true out=HD_to_seq_filtered_beagle5
tabix -f HD_to_seq_filtered_beagle5.vcf.gz
```

```
beagle gt=study_hd.vcf.gz ref=ref_nonfiltered_phased.vcf.gz chrom=13 impute=true gp=true out=HD_to_seq_nonfiltered_beagle5
tabix -f HD_to_seq_nonfiltered_beagle5.vcf.gz
```

## 8. Imputation using minimac3

The imputation process for using minimac3 is rather similar. It is more efficient than Beagle 5 but slightly slower. It takes ~10 mins to impute to filtered sequence reference and ~15 mins to impute to unfiltered sequence reference. I have already done the process, so you can just copy the outputs from the project folder to the current imputation folder. 

```
#$minimac3 --refHaps ref_filtered_phased.vcf.gz --haps study_hd.vcf.gz --prefix HD_to_seq_filtered_minimac3
#tabix -f HD_to_seq_filtered_minimac3.dose.vcf.gz
#$minimac3 --refHaps ref_nonfiltered_phased.vcf.gz --haps study_hd.vcf.gz --prefix HD_to_seq_nonfiltered_minimac3
#tabix -f HD_to_seq_nonfiltered_minimac3.dose.vcf.gz
```

```
cp $maindir/HD_to_seq_filtered_minimac3.* $HOME/imputation-workshop/imputation
cp $maindir/HD_to_seq_nonfiltered_minimac3.* $HOME/imputation-workshop/imputation
```

## 9. Calculate the genotype concordance using vcf-compare (from VCFtools)

In this tutorial, I am gonna show you two parameters: genotype concordance and allelic/dosage R-square.

To compare two vcfs and have an idea of genotype concordance, there is a sub-function from vcftools: vcf-compare. so just pop in vcf-compare VCF1 VCF2 > output. You will have an output file.

In the previous session, we have four imputation output using both Beagle 5 and minimac3 to impute to filtered and unfiltered sequence reference. So four concordance file will be generated as below:

```
vcf-compare study_filtered.vcf.gz HD_to_seq_filtered_beagle5.vcf.gz > concordance_beagle5_filtered
less concordance_beagle5_filtered
```
![](https://github.com/GenomicsAotearoa/Imputation-workshop/blob/master/Tutorial/img/Screen%20Shot%202020-09-16%20at%2014.50.01.png?raw=true)

```
vcf-compare study_filtered.vcf.gz HD_to_seq_nonfiltered_beagle5.vcf.gz > concordance_beagle5_nonfiltered`
less concordance_beagle5_nonfiltered
```
![](https://github.com/GenomicsAotearoa/Imputation-workshop/blob/master/Tutorial/img/Screen%20Shot%202020-09-16%20at%2014.51.21.png?raw=true)

```
vcf-compare study_filtered.vcf.gz HD_to_seq_filtered_minimac3.dose.vcf.gz > concordance_minimac3_filtered
less concordance_minimac3_filtered
```
![](https://github.com/GenomicsAotearoa/Imputation-workshop/blob/master/Tutorial/img/Screen%20Shot%202020-09-16%20at%2014.52.43.png?raw=true)

```
vcf-compare study_filtered.vcf.gz HD_to_seq_nonfiltered_minimac3.dose.vcf.gz > concordance_minimac3_nonfiltered
less concordance_minimac3_nonfiltered
```
![](https://github.com/GenomicsAotearoa/Imputation-workshop/blob/master/Tutorial/img/Screen%20Shot%202020-09-16%20at%2014.53.47.png?raw=true)

## 10. Evaluate the performance of imputation: allelic/dosage R-square

To calculate the dosage R-square, beagle 5 does not provide a seperate file. You may need a bit code to extract the information:

```
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%DR2\t%AF\t%IMP\n' HD_to_seq_filtered_beagle5.vcf.gz > HD_to_seq_filtered_beagle5.r2
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%DR2\t%AF\t%IMP\n' HD_to_seq_nonfiltered_beagle5.vcf.gz > HD_to_seq_nonfiltered_beagle5.r2
```

The columns of the file we generated are chromosome, position, SNP name, reference allele, alternative allele, quality, filter, **dosage r-square,** allele frequency, whether it is imputed. It is a thoughtful enough file that provides us all the information, the only additional part we may need to do is calculate the minor allele frequency from allele frequency.

```
head HD_to_seq_filtered_beagle5.r2
```
![](https://github.com/GenomicsAotearoa/Imputation-workshop/blob/master/Tutorial/img/Screen%20Shot%202020-09-16%20at%2014.56.35.png?raw=true)

minimac3 generate an info file after it finishes imputing. It is pretty thoughtful that it provides us minor allele frequency directly. The troubling part is that we have to extract the position from the SNP column for future comparison. 

```
head HD_to_seq_filtered_minimac3.info
```

![](https://github.com/GenomicsAotearoa/Imputation-workshop/blob/master/Tutorial/img/Screen%20Shot%202020-09-16%20at%2015.02.07.png?raw=true)

So we got four output here. and we are gonna pop them in R to have a look: 

`R`

The first step is to read in all our output files in R

```
filteredBG5 <- read.table("HD_to_seq_filtered_beagle5.r2")
nonfilteredBG5 <- read.table("HD_to_seq_nonfiltered_beagle5.r2")
filteredminimac3 <- read.table("HD_to_seq_filtered_minimac3.info", header=T)
nonfilteredminimac3 <- read.table("HD_to_seq_nonfiltered_minimac3.info", header=T)
```

Then we need to extract all the positions. This step is a bit redundant for beagle outputs but really helpful for the minimac3 output. The function we are gonna use is `substr`, it tells R to just extract the string from the 4th digit to the 11th digit. 
 
```
filteredBG5$Pos <- filteredBG5$V2
nonfilteredBG5$Pos <- nonfilteredBG5$V2
filteredminimac3$Pos <- substr(filteredminimac3$SNP, 4, 11)
nonfilteredminimac3$Pos <- substr(nonfilteredminimac3$SNP, 4, 11)
```
The next step is to extract the R-square for beagle 5. Usually, it shouldn't be a problem, you get the number in that column directly. However, in this session, we used the unfiltered reference, which contains the multi-allelic positions. In this case, Beagle will give you multiple possible solutions for those multi-allelic positions. In this case, we just take the first solution to make things easier. Here you may see a warning message mentioned `NA` generated. Don't worry about that.  

```r
filteredBG5$DR2_filtered_BG5 <- filteredBG5$V8
nonfilteredBG5$DR2_nonfiltered_BG5 <- as.numeric(substr(nonfilteredBG5$V8, 1, 4))
filteredminimac3$Rsq_filtered_minimac3 <- filteredminimac3$Rsq
nonfilteredminimac3$Rsq_nonfiltered_minimac3 <- nonfilteredminimac3$Rsq
```
Now let's merge both the output files from Beagle and Minimac3, then final merge them into a file called `finalmerge`

```r
mergedBeagle <- merge(filteredBG5, nonfilteredBG5, by.x="Pos", by.y="Pos", all=FALSE)
mergedMinimac3 <- merge(filteredminimac3, nonfilteredminimac3, by.x="Pos", by.y="Pos", all=FALSE)
finalmerge <- merge(mergedBeagle, mergedMinimac3, by.x="Pos", by.y="Pos", all=FALSE)
```

Let's have a look at the summary

```r
summary(finalmerge$DR2_filtered_BG5)
summary(finalmerge$DR2_nonfiltered_BG5)
summary(finalmerge$Rsq_filtered_minimac3)
summary(finalmerge$Rsq_nonfiltered_minimac3)
```

In both cases of beagle 5 and minimac3, using unfiltered reference gave us poorer performance compared to the filtered ones. Minimac3 gave slightly higher allelic square compare to beagle5. It is not always the case since we are only using 5MB here. And also the performance depends on a lot of parameters. As I mentioned, Beagle is fast but computationally demanding. Minimac 3 is slower but very efficient. There are of course other software for you to choose. Which software to use, what parameters for QC, questions such as those I may not have an answer, you have to figure it out by doing experiments. 

As I also mentioned allelic/dosage-r square is a good parameter for evaluating the performance. Here you can see the relationship between MAF and allelic/dosage-r square. 

![](https://github.com/GenomicsAotearoa/Imputation-workshop/blob/master/Tutorial/img/Screen%20Shot%202020-09-16%20at%2015.13.24.png?raw=true)

```
pdf()
plot(finalmerge$MAF.x,finalmerge$DR2_filtered_BG5, pch=4)
dev.off()
```

![](https://github.com/GenomicsAotearoa/Imputation-workshop/blob/master/Tutorial/img/image003.png?raw=true)

```
pdf()
plot(finalmerge$MAF.x,finalmerge$Rsq_filtered_minimac3, pch=4)
dev.off()
```

![](https://github.com/GenomicsAotearoa/Imputation-workshop/blob/master/Tutorial/img/image002.png?raw=true)

And also, we can also have a look at the correlation between the allelic/dosage-r square from beagle 5 and minimac3. 

```
pdf()
plot(finalmerge$DR2_filtered_BG5,finalmerge$Rsq_filtered_minimac3, pch=4)
dev.off()
```

![](https://github.com/GenomicsAotearoa/Imputation-workshop/blob/master/Tutorial/img/image001.png?raw=true)
