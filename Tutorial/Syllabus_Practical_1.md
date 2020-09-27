# OVERVIEW
In imputation studies, there are two datasets involved: a reference and a target. The reference dataset is a set of individuals that have genotype data at the desired density, while the target data includes individuals that you wish to infer the higher density genotypes. A crucial aspect to any imputation study is identifying the appropriate reference set of individuals to represent the target individuals. However, scenario where only a limited subset of individuals can be genotyped (e.g. due to cost) at the desired density, choosing the right individuals to serve as a reference set become critically important.

# AIMS:
1)	Assess the impact that using different individuals in a reference set has on imputation accuracy.

# BACKGROUND
What has already been done to the data before this point.
Including running LP Choose and refer to code you will add in.

# TOOLS:
R: Packages data.table and reshape2 – offer quicker manipulations of large data frames

# INPUT FILES:
- infoFile.csv: Breed information on each animal as well as Principal Components from full genome
- TrueGenotypes.vcf: True Genotypes on each animal to compute accuracy statistics
- ImputedData_xxx.vcf: Imputed genotypes based on reference population xxx with options
   - xxx1
   - xxx2
   - xxx3
   
# PROCEDURES:
## 1) Load modules
For this workshop we only need to use the R module because the imputation has already been completed.
```
cd $HOME/imputation-workshop

module load R
```

## 2) Generate PCA Plots in R
Open R 
```
R
```

Generate link to the workshop folder
```
maindir = "/nesi/nobackup/nesi02659/SEP28/practical2/"
```

Read in the info file and take a look at it
```
infoFile = paste(maindir,"infoFile.csv",sep="")
head(infoFile)
```

A PCA plot can give us an idea of how closely related the reference and target animals are. I have already generated the PCs using the full genome (rather than the single chromosome used in this example) and this was done using the "eigen" command run on the genomic relationship matrix. This code will be made available on the github site if you are interested in working through it independently.

The "PlotColour" column has been set for different colours based on breed, with WGS animals a different colour, based on either "NZ" or "Other" as the country that provided the WGS.
```
par(bg = "gray") # set a gray background so you can see the colours more clearly
plot(infoFile$PC1,infoFile$PC2,pch=19,col=infoFile$PlotColour)
```

Now that you have the PCA plot you can look at the relationships between those animals that have WGS and will be used as the reference set of animals, and the target population.

Some other things you can try by changing the "PlotColour" column in the above code"
- How do Australian WGS animals compare with Other WGS animals?
- How representative are the samples that were selected by LPChoose (LPChoose column)? These were selected based on the WGS animals already being in the reference population.

## 3) Impact of Different Imputation Sets


STEP 6: Assess the impact of using different reference sets for imputation

Background: Here, we’ll test the impact of using different reference sets for imputation of individuals with 50k genotypes to HD genotypes. We have created a series of files containing different sets of reference individuals:

LPChoose_HD.ans
Random_HD.ans
UNREL_WGS_IDs.txt
UNREL_AUS_WGS_IDs.txt
ALL_WGS_IDs.txt
NZ_WGS_IDs.txt
NZ_AUS_WGS_IDs.txt
AUS_WGS_IDs.txt

```
These files can be used to subset our HD genotype file HD_SNPs_26_All_Ans.vcf.gz, e.g. in bcftools:
bcftools view -Oz -S UNREL_WGS_IDs.txt HD_SNPs_26_All_Ans.vcf.gz > HD_SNPs_26_Unrelated_WGS_Ans.vcf.gz

bcftools view -Oz -S NZ_WGS_IDs.txt HD_SNPs_26_All_Ans.vcf.gz > HD_SNPs_26_NZ_WGS_Ans.vcf.gz

We will also want to make sure our 50k genotypes don’t contain any reference individuals:
bcftools view -Oz -S Target_HD_50K_Overlap.ans HD_50K_Overlap_SNPs_26_All_Ans.vcf.gz > HD_50K_Overlap_SNPs_26_Target_Ans.vcf.gz

To assess imputation accuracy, we have masked 10% of the 50k genotypes:
HD_50K_Overlap_Masked_SNPs.pos
HD_50K_Overlap_Kept_SNPs.pos

##50k SNP set for Imputation##
bcftools view -Oz -T HD_50K_Overlap_Kept_SNPs.pos HD_50K_Overlap_SNPs_26_Target_Ans.vcf.gz > HD_50K_Overlap_SNPs_26_Target_Ans_kept.vcf.gz

##50k SNP masked set of genotypes – true genotypes##
bcftools view -Oz -T HD_50K_Overlap_Masked_SNPs.pos HD_50K_Overlap_SNPs_26_Target_Ans.vcf.gz > HD_50K_Overlap_SNPs_26_Target_Ans_masked.vcf.gz

gunzip HD_50K_Overlap_SNPs_26_Target_Ans_masked.vcf.gz

sed -i 's@0/0:.:.:.:.@'0'@g' HD_50K_Overlap_SNPs_26_Target_Ans_masked.vcf
sed -i 's@0/1:.:.:.:.@'1'@g' HD_50K_Overlap_SNPs_26_Target_Ans_masked.vcf
sed -i 's@1/0:.:.:.:.@'1'@g' HD_50K_Overlap_SNPs_26_Target_Ans_masked.vcf
sed -i 's@1/1:.:.:.:.@'2'@g' HD_50K_Overlap_SNPs_26_Target_Ans_masked.vcf
sed -i 's@./.:.:.:.:.@'NA'@g' HD_50K_Overlap_SNPs_26_Target_Ans_masked.vcf

From here, we can impute using various reference sets, e.g.:
##Phase reference sets##
beagle gt=HD_SNPs_26_Unrelated_WGS_Ans.vcf.gz nthreads=16 ne=500 out=HD_SNPs_26_Unrelated_WGS_Ans_phased

beagle gt=HD_SNPs_26_NZ_WGS_Ans.vcf.gz nthreads=16 ne=500 out=HD_SNPs_26_NZ_WGS_Ans_phased

##Impute target using reference##
beagle gt=HD_50K_Overlap_SNPs_26_Target_Ans_kept.vcf.gz ref=HD_SNPs_26_Unrelated_WGS_Ans_phased.vcf.gz nthreads=16 ne=500 out=HD_50K_Overlap_SNPs_26_Target_Ans_Imputed_with_Unrelated_WGS_Ans

beagle gt=HD_50K_Overlap_SNPs_26_Target_Ans_kept.vcf.gz ref=HD_SNPs_26_NZ_WGS_Ans_phased.vcf.gz nthreads=16 ne=500 out=HD_50K_Overlap_SNPs_26_Target_Ans_Imputed_with_NZ_WGS_Ans

To assess imputation accuracy, we’ll trim down the output to only include the masked SNPs:
bcftools index HD_50K_Overlap_SNPs_26_Target_Ans_Imputed_with_Unrelated_WGS_Ans.vcf.gz
bcftools view -Oz -T HD_50K_Overlap_Masked_SNPs.pos HD_50K_Overlap_SNPs_26_Target_Ans_Imputed_with_Unrelated_WGS_Ans.vcf.gz > HD_50K_Overlap_SNPs_26_Target_Ans_Imputed_with_Unrelated_WGS_Ans_masked_SNPs.vcf.gz

bcftools index HD_50K_Overlap_SNPs_26_Target_Ans_Imputed_with_NZ_WGS_Ans.vcf.gz
bcftools view -Oz -T HD_50K_Overlap_Masked_SNPs.pos HD_50K_Overlap_SNPs_26_Target_Ans_Imputed_with_NZ_WGS_Ans.vcf.gz > HD_50K_Overlap_SNPs_26_Target_Ans_Imputed_with_NZ_WGS_Ans_masked_SNPs.vcf.gz

```

##ASSESS IMPUTATION ACCURACY FOR THE DIFFERENT SETS - CHANGE "imputed_file"###

```{R}
require(data.table)
library(reshape2)

imputed_file1="HD_50K_Overlap_SNPs_26_Target_Ans_Imputed_with_Unrelated_WGS_Ans_masked_SNPs.vcf.gz"
imputed_file2="HD_50K_Overlap_SNPs_26_Target_Ans_Imputed_with_NZ_WGS_Ans_masked_SNPs.vcf.gz"
true=fread("HD_50K_Overlap_SNPs_26_Target_Ans_masked.vcf",skip=5765)
true2=true[,-c(1:9)]
true2=as.matrix(true2)

imputed1=fread(imputed_file1)
imputed2=fread(imputed_file2)

imputeda=imputed1[,-c(1:9)]
imputeda=as.data.frame(unlist(imputeda))
colnames(imputeda)="var"
imputeda <- colsplit(imputeda$var, ":", c("left", "right"))
a=matrix(imputeda$left,nrow(imputed1),ncol(imputed1)-9)
a[a=="1|1"]<-2
a[a=="1|0"]<-1
a[a=="0|1"]<-1
a[a=="0|0"]<-0
a=as.data.frame(a)
a2=matrix(imputeda$right,nrow(imputed1),ncol(imputed1)-9)

imputedb=imputed2[,-c(1:9)]
imputedb=as.data.frame(unlist(imputedb))
colnames(imputedb)="var"
imputedb <- colsplit(imputedb$var, ":", c("left", "right"))
b=matrix(imputedb$left,nrow(imputed2),ncol(imputed2)-9)
b[b=="1|1"]<-2
b[b=="1|0"]<-1
b[b=="0|1"]<-1
b[b=="0|0"]<-0
b=as.data.frame(b)
b2=matrix(imputedb$right,nrow(imputed2),ncol(imputed2)-9)


##Concordance##
al=a==true2
al2=al
al2[which(is.na(true2))] = NA
acca=colMeans(al2,na.rm=T) ##Individual Concordance
sncona=rowMeans(al2,na.rm=T)##SNP Concordance
summary(sncona)

bl=b==true2
bl2=bl
bl2[which(is.na(true2))] = NA
accb=colMeans(bl2,na.rm=T) ##Individual Concordance
snconb=rowMeans(bl2,na.rm=T)##SNP Concordance
summary(snconb)

##Dosage R2 from Beagle##
dra=sapply(strsplit(as.character(imputed1$INFO),";"),"[",1)
dra=sapply(strsplit(dra,"="),"[",2)
dra=as.numeric(dra)
summary(dra)

drb=sapply(strsplit(as.character(imputed2$INFO),";"),"[",1)
drb=sapply(strsplit(drb,"="),"[",2)
drb=as.numeric(drb)
summary(drb)

##SNP accuracy from absolute calls##
asncora=cor(as.numeric(a[1,]),as.numeric(true2[1,]),use="complete.obs")
for(i in 2:nrow(a)){
	asncora=c(asncora,cor(as.numeric(a[i,]),as.numeric(true2[i,]),use="complete.obs"))
}
summary(asncora)

asncorb=cor(as.numeric(b[1,]),as.numeric(true2[1,]),use="complete.obs")
for(i in 2:nrow(b)){
	asncorb=c(asncorb,cor(as.numeric(b[i,]),as.numeric(true2[i,]),use="complete.obs"))
}
summary(asncorb)

##SNP accuracy from dosage calls##
dsncora=cor(as.numeric(a2[1,]),as.numeric(true2[1,]),use="complete.obs")
for(i in 2:nrow(a2)){
	dsncora=c(dsncora,cor(as.numeric(a2[i,]),as.numeric(true2[i,]),use="complete.obs"))
}
summary(dsncora)

dsncorb=cor(as.numeric(b2[1,]),as.numeric(true2[1,]),use="complete.obs")
for(i in 2:nrow(b2)){
	dsncorb=c(dsncorb,cor(as.numeric(b2[i,]),as.numeric(true2[i,]),use="complete.obs"))
}
summary(dsncorb)

##Individual accuracy from absolute calls##
aancora=cor(as.numeric(unlist(a[,1])),as.numeric(unlist(true2[,1])),use="complete.obs")
for(i in 2:ncol(a)){
	aancora=c(aancora,cor(as.numeric(unlist(a[,i])),as.numeric(unlist(true2[,i])),use="complete.obs"))
}
summary(aancora)

aancorb=cor(as.numeric(unlist(b[,1])),as.numeric(unlist(true2[,1])),use="complete.obs")
for(i in 2:ncol(b)){
	aancorb=c(aancorb,cor(as.numeric(unlist(b[,i])),as.numeric(unlist(true2[,i])),use="complete.obs"))
}
summary(aancorb)

##Individual accuracy from dosage calls##
dancora=cor(as.numeric(unlist(a2[,1])),as.numeric(unlist(true2[,1])),use="complete.obs")
for(i in 2:ncol(a2)){
	dancora=c(dancora,cor(as.numeric(unlist(a2[,i])),as.numeric(unlist(true2[,i])),use="complete.obs"))
}
summary(dancora)

dancorb=cor(as.numeric(unlist(b2[,1])),as.numeric(unlist(true2[,1])),use="complete.obs")
for(i in 2:ncol(a2)){
	dancorb=c(dancorb,cor(as.numeric(unlist(b2[,i])),as.numeric(unlist(true2[,i])),use="complete.obs"))
}
summary(dancorb)

##Compare reference sets##
t.test(dancora, dancorb, paired = TRUE, alternative = "two.sided")

```
