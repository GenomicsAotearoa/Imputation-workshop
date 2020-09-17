# Aim
1. Go through the pipeline of phasing and imputing high density genotype to sequence level
2. Understand the importance of quality control 
3. Know how to evaluate the imputation performance

# Tools we need
BCFtools: basic bioinformatics software, in this tutorial, we use it for creating subsets and quality control ([http://samtools.github.io/bcftools/bcftools.html](http://samtools.github.io/bcftools/bcftools.html))

VCFtools: basic bioinformatics software, in this tutorial, we use it for compare two vcf files and evaluate the concordance rate ([https://vcftools.github.io/index.html](https://vcftools.github.io/index.html))

R: basic statistics software

Beagle 5.0: software for phasing and imputation ([https://faculty.washington.edu/browning/beagle/beagle.html](https://faculty.washington.edu/browning/beagle/beagle.html))

Minimac3: software for imputation ([https://genome.sph.umich.edu/wiki/Minimac3](https://genome.sph.umich.edu/wiki/Minimac3))

# Input files
In this tutorial, we used 5MB segment (30-35MB) of chr13 from the 1000 Genome phase 3, which is publicly available and can be downloaded from https://mathgen.stats.ox.ac.uk/impute/impute_v2.html

# Procedures
## 1. Load the module that required on NeSI
To find modules which are already available on NeSI, use `module spider #module_name`

To load modules and start using, use `module load #module_name`

```bash
module load BCFtools
module load VCFtools
module load R
beagle5= #(TO BE TRANSFERED TO TRAINING PROJECT FOLDER)
minimac3= #(TO BE TRANSFERED TO TRAINING PROJECT FOLDER)
```

## 2. Copy the files into the home directory
Define two directories: workshop directory and home directory. In this workshop, the analysis will be conducted in the home directory

```bash
maindir=/scale_wlg_nobackup/filesets/nobackup/nesi00187/yuwan0/workshop
home=/home/ywang17
```

The main input file (seqvcf) is extracted (from 30MB to 35MB) on chr13 from 1000 human genome data. I selected some of the reliable SNPs to generate a HD genotype dataset(hdvcf). 

```
seqvcf=$maindir/nonfilter_seq_5MB.vcf.gz
hdvcf=$maindir/hd_5MB.vcf.gz
```

Now what you need to do is to cp both genotype and sequence data to your own home directory. In addition, we have to software that we need to use which are not avaiable on NeSI. Please also copy these two files in your own working directory.

```bash
cp $seqvcf $home
cp $hdvcf $home
cp $beagle5 $home
cp $minimac3 $home
beagle5=/home/ywang17/beagle.12Jul19.0df_bg5.jar
minimac3=/home/ywang17/Minimac3
```

The genotype and sequence files use "vcf.gz" format. We can not open it directly. To check how the genotype looks like, you need to use: zless -S. -S is to make the file well formated. 

```bash
zless -S $seqvcf
zless -S $hdvcf
```

## 3. Explore the input file

To have a basic idea of the genotype, BCFtools have a very convenient function: `stats`. By checking the original sequence file's information. the code you need is as below. -s is a common tag to show "samples". Samples to include or "-" to apply all variants. Via adding this, we also generate the statistics for each individual. The output will be named: original.stats

`bcftools stats -s "-" nonfilter_seq_5MB.vcf.gz > original.stats`

Now let's have a look at the output:

`less original.stats`

So in the first part you can see the basic statistics of the sequence file. We have 2504 samples in total, with 149,854 variants. 144,069 of them are SNPs, 5777 indels, 98 other, 694 multi-allelic sites with 467 multi-allelic SNPs.  

## 4. Quality control process

There are a lot of parameters you may take into consideration in your dataset, such as non-variant, singletons, multi-allelic positions, map quality, mendelian error, minor allele frequency etc. The parameter settings vary depends highly on data, and also for the purpose of data. For example, if the data is for GWAS analysis, you probably don't want to throw away too much variants since a lot of the causals are rare alleles. If your data is used for prediction, then you probably don't want to include too much crap in the data, which is not beneficial for the next step model construction.

In this tutorial, I am only gonna consider non-variant, singleton and multi-allelic variants since majority of the imputation software can not handle them anyway. Other filtering process can also be done using bcftools, just pop in the website and go to filtering session.

```bash
bcftools view -O z -o nosingleton.vcf.gz -i 'AC>1' nonfilter_seq_5MB.vcf.gz
tabix -f nosingleton.vcf.gz
bcftools stats -s "-" nosingleton.vcf.gz > step1.stats
```

`less step1.stats`

```bash
bcftools view -O z -o nosingleton_2alleles.vcf.gz --max-alleles 2 nosingleton.vcf.gz
tabix -f nosingleton_2alleles.vcf.gz
bcftools stats -s "-" nosingleton_2alleles.vcf.gz > step2.stats
```

```bash
less step2.stats
```

## 5. Generate the data for the imputation process
We are gonna just use this data for imputation. So there are two populations involved: reference population and study population. I decide to use the first 1000 samples as my reference, and the rest will be set as my study population. Using bcftools to extract the sample ID and awk to generate two population ID files as follow:

```bash
bcftools query -l nonfilter_seq_5MB.vcf.gz > seq.ID
awk 'NR>0&&NR<=1000' seq.ID > ref.ID
awk 'NR>1000&&NR<=2000' seq.ID > study.ID
```

Now we create a new directory for our imputation :

```bash
mkdir -p imputation
cd imputation
```

bcftools is very handy of extracting the sample from the whole dataset via using -S. We will generate two files for our study samples: one from HD genotype as study. And another one from filtered sequence data for future validation. After extraction, we use function tabix to generic indexer for TAB-delimited genome position files.

```bash
bcftools view -O z -o study_hd.vcf.gz -S $home/study.ID $home/hd_5MB.vcf.gz
tabix -f study_hd.vcf.gz
bcftools view -O z -o study_filtered.vcf.gz -S $home/study.ID $home/nosingleton_2alleles.vcf.gz #for validation
tabix -f study_filtered.vcf.gz
```

Similarly we will extract the samples for our reference samples. Same functions will be used here. We will generate two files for the reference samples, one is from the original unfiltered sequence, and another one from filtered sequence. Again, we will tabix these two files. 

```bash
bcftools view -O z -o ref_nonfiltered.vcf.gz -S $home/ref.ID $home/nonfilter_seq_5MB.vcf.gz
tabix -f ref_nonfiltered.vcf.gz
bcftools view -O z -o ref_filtered.vcf.gz -S $home/ref.ID $home/nosingleton_2alleles.vcf.gz
tabix -f ref_filtered.vcf.gz
```

## 6. Phasing using Beagle 5 

For imputation, no matter which software are you using, phasing is compulsory for the reference population. For some software like minimac3, both study and reference population have to be phased.  If you simply pop the unphased reference population into your imputation software, the software will immediately give you an error message. 

The data I downloaded already finished phasing that you can see in the dataset, all the genotypes are phased (eg: "0|1", not "0/1"). Besides phasing needs a lot of computation resources and certain amount of time. So I will just show you the code for the phasing process.    

```bash
#java -jar $beagle5 gt=ref_filtered.vcf.gz chrom=13 out=ref_filtered_phased
#tabix -f ref_filtered_phased.vcf.gz
```

```bash
#java -jar $beagle5 gt=ref_nonfiltered.vcf.gz chrom=13 out=ref_nonfiltered_phased
#tabix -f ref_nonfiltered_phased.vcf.gz
```

## 7. Imputation using Beagle 5

In this tutorial, I will show you the imputation using two software: Beagle 5 and minimac3. Both software are very stable, reliable, easy-to-use, free, and pretty popular. Beagle 5 is computationally demanding but can give you the accurate results very fast. Minimac is computationally efficient, but a bit slower. In addition, to prove that quality control is an important procedure, both filtered reference and unfiltered sequence reference will be used.

Beagle has been evolved from version 3.0 to current 5.1. The parameters for running the software have been reduced a lot. There are several important parameters that can influence imputation performance such as effective population size (Ne), window size, etc. Check the following paper: Improving Imputation Quality in BEAGLE for Crop and Livestock Data [https://www.g3journal.org/content/10/1/177](https://www.g3journal.org/content/10/1/177)

```bash
java -jar $beagle5 gt=study_hd.vcf.gz ref=ref_filtered_phased.vcf.gz chrom=13 impute=true gp=true out=HD_to_seq_filtered_beagle5
#Imputation time: 54 seconds
#Total time: 57 seconds
```

```bash
java -jar $beagle5 gt=study_hd.vcf.gz ref=ref_nonfiltered_phased.vcf.gz chrom=13 impute=true gp=true out=HD_to_seq_nonfiltered_beagle5
#Imputation time: 1 minute 21 seconds
#Total time: 1 minute 25 seconds
```

```bash
tabix -f HD_to_seq_filtered_beagle5.vcf.gz
tabix -f HD_to_seq_nonfiltered_beagle5.vcf.gz
```

## 8. Imputation using minimac3

The imputation process for using minimac3 is rather similar. It takes ~10 mins to impute to filtered sequence reference and ~15 mins to impute to unfiltered sequence reference. I have already done the process, so you can just copy the outputs from the project folder to current imputation folder. 

```bash
#$minimac3 --refHaps ref_filtered_phased.vcf.gz --haps study_hd.vcf.gz --prefix HD_to_seq_filtered_minimac3
```

```bash
#tabix -f HD_to_seq_filtered_minimac3.dose.vcf.gz
```

```bash
#$minimac3 --refHaps ref_nonfiltered_phased.vcf.gz --haps study_hd.vcf.gz --prefix HD_to_seq_nonfiltered_minimac3
```

```bash
#tabix -f HD_to_seq_nonfiltered_minimac3.dose.vcf.gz
```

```bash
cp $maindir/HD_to_seq_filtered_minimac3.* $home/imputation
cp $maindir/HD_to_seq_nonfiltered_minimac3.* $home/imputation
```

## 9. Calculate the genotype concordance using vcf-compare (from VCFtools)

In this tutorial, I am gonna show you two parameter: genotype concordance and allelic/dosage R-square.

To compare two vcfs and have an idea of genotype concordance, there is a sub function from vcftools: vcf-compare. so just pop in vcf-compare VCF1 VCF2 > output. You will have a output.

In the previous session, we have four imputation output using both Beagle 5 and minimac to impute to filtered and unfiltered sequence reference. So four concordance file will be generated as below:

```bash
vcf-compare study_filtered.vcf.gz HD_to_seq_filtered_beagle5.vcf.gz > concordance_beagle5_filtered
less concordance_beagle5_filtered
```


```bash
vcf-compare study_filtered.vcf.gz HD_to_seq_nonfiltered_beagle5.vcf.gz > concordance_beagle5_nonfiltered
less concordance_beagle5_nonfiltered
```


```bash
vcf-compare study_filtered.vcf.gz HD_to_seq_filtered_minimac3.dose.vcf.gz > concordance_minimac3_filtered
less concordance_minimac3_filtered
```


```bash
vcf-compare study_filtered.vcf.gz HD_to_seq_nonfiltered_minimac3.dose.vcf.gz > concordance_minimac3_nonfiltered
less concordance_minimac3_nonfiltered
```



## 10. Evaluate the performance of imputation: allelic/dosage R-square

To calculate the dosage R-square, beagle 5 does not provide a seperate file. You may need a bit code to extract the information:

```bash
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%DR2\t%AF\t%IMP\n' HD_to_seq_filtered_beagle5.vcf.gz > HD_to_seq_filtered_beagle5.r2
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%DR2\t%AF\t%IMP\n' HD_to_seq_nonfiltered_beagle5.vcf.gz > HD_to_seq_nonfiltered_beagle5.r2
```

The columns of the file we generated are chromosome, position, SNP name, reference allele, alternative allele, quality, filter, **dosage r-square,** allele frequency, whether it is imputed. It is a thoughtful enough file that provides us all the information, the only additional part we may need to do is calculate the minor allele frequency from allele frequency.

```bash
head HD_to_seq_filtered_beagle5.r2
```



```bash
head HD_to_seq_filtered_minimac3.info
```

minimac generate an info file after it finishes imputing. It is pretty thoughtful that it provides us minor allele frequency directly. The trouble part is that we have to extract the position from SNP column for future compare. 


So we got four output here. and we are gonna pop them in R to have a look: 

```r
filteredBG5 <- read.table("HD_to_seq_filtered_beagle5.r2")
nonfilteredBG5 <- read.table("HD_to_seq_nonfiltered_beagle5.r2")
filteredminimac3 <- read.table("HD_to_seq_filtered_minimac3.info", header=T)
nonfilteredminimac3 <- read.table("HD_to_seq_nonfiltered_minimac3.info", header=T)
```

```r
filteredBG5$Pos <- filteredBG5$V2
nonfilteredBG5$Pos <- nonfilteredBG5$V2
filteredminimac3$Pos <- substr(filteredminimac3$SNP, 4, 11)
nonfilteredminimac3$Pos <- substr(nonfilteredminimac3$SNP, 4, 11)
```

```r
filteredBG5$DR2_filtered_BG5 <- filteredBG5$V8
nonfilteredBG5$DR2_nonfiltered_BG5 <- as.numeric(substr(nonfilteredBG5$V8, 1, 4))
filteredminimac3$Rsq_filtered_minimac3 <- filteredminimac3$Rsq
nonfilteredminimac3$Rsq_nonfiltered_minimac3 <- nonfilteredminimac3$Rsq
```

```r
mergedBeagle <- merge(filteredBG5, nonfilteredBG5, by.x="Pos", by.y="Pos", all=FALSE)
mergedMinimac3 <- merge(filteredminimac3, nonfilteredminimac3, by.x="Pos", by.y="Pos", all=FALSE)
finalmerge <- merge(mergedBeagle, mergedMinimac3, by.x="Pos", by.y="Pos", all=FALSE)
```

```r
summary(finalmerge$DR2_filtered_BG5)
summary(finalmerge$DR2_nonfiltered_BG5)
summary(finalmerge$Rsq_filtered_minimac3)
summary(finalmerge$Rsq_nonfiltered_minimac3)
```


`pdf()`
`plot(finalmerge$MAF.x,finalmerge$DR2_filtered_BG5, pch=4)`
`dev.off()`

`pdf()`
`plot(finalmerge$MAF.x,finalmerge$Rsq_filtered_minimac3, pch=4)`
`dev.off()`

