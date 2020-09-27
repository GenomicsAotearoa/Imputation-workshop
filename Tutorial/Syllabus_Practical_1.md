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
 '/nesi/nobackup/nesi02659/SEP28/practical2/infoFile.csv'
```

A PCA plot can give us an idea of how closely related the reference and target animals are. I have already generated the PCs using the full genome (rather than the single chromosome used in this example) and this was done using the "eigen" command run on the genomic relationship matrix. This code will be made available on the github site if you are interested in working through it independently.

The "PlotColour" column has been set for different colours based on breed, with WGS animals a different colour, based on either "NZ" or "Other" as the country that provided the WGS.
```
par(bg = "gray") # set a gray background so you can see the colours more clearly
plot(infoFile$PC1,infoFile$PC2,pch=19,col=infoFile$PlotColour)
```

Now that you have the PCA plot you can look at the relationships between those animals that have WGS and will be used as the reference set of animals, and the target population.

#### Extra Exploration of Data
Some other things you can try by changing the "PlotColour" column in the above code"
- How do Australian WGS animals compare with Other WGS animals?
- How representative are the samples that were selected by LPChoose (LPChoose column)? These were selected based on the WGS animals already being in the reference population.

## 3) Impact of Different Imputation Sets
Here, we’ll test the impact of using different reference sets for imputation of individuals with 50k genotypes to HD genotypes. We have run BEAGLE 5.1 a number of times using a variety of reference populations and provided the output files for you to evaluate imputation accuracy. These files contain the imputed genotypes for the target animals (i.e. animals with 50k genotypes) for a set of SNPs that have been masked to test imputation accuracy.

Reference Population:
- All WGS animals: ImputedData_ALL_WGS.vcf.gz
- NZ WGS animals: ImputedData_NZ_WGS.vcf.gz
- NZ + AUS WGS animals: ImputedData_NZ_AUS_WGS.vcf.gz
- AUS WGS animals: ImputedData_AUS_WGS.vcf.gz
- Other WGS animals: ImputedData_Other_WGS.vcf.gz
- Other + AUS animals: ImputedData_Other_AUS_WGS.vcf.gz
- LPChoose + WGS animals: ImputedData_LPChoose_WGS.vcf.gz
- Random + WGS animals: ImputedData_Random_WGS.vcf.gz

First we will call 2 packages that will help read and format large datasets
```
require(data.table)
library(reshape2)
```

Then we will load some functions to help us process the data
```
getGenotypeCalls = function(imputedVCF){
   imputed= imputedVCF[,-c(1:9)]
   imputed=as.data.frame(unlist(imputed))
   colnames(imputed)="var"
   imputed <- colsplit(imputed$var, ":", c("genotype", "dosage"))
   a=matrix(imputed$genotype,nrow(imputedVCF))
   a[a=="1|1"]<-2
   a[a=="1|0"]<-1
   a[a=="0|1"]<-1
   a[a=="0|0"]<-0
   a=as.data.frame(a)
   colnames(a) = colnames(imputedVCF)[10:ncol(imputedVCF)]
   rownames(a) = paste("C", imputedVCF$"#CHROM","_P", imputedVCF$POS,sep="")
   a
}

getDosages = function(imputedVCF){
   imputed= imputedVCF[,-c(1:9)]
   imputed=as.data.frame(unlist(imputed))
   colnames(imputed)="var"
   imputed <- colsplit(imputed$var, ":", c("genotype", "dosage"))
   a=as.data.frame(matrix(imputed$dosage,nrow(imputedVCF)))
   colnames(a) = colnames(imputedVCF)[10:ncol(imputedVCF)]
   rownames(a) = paste("C", imputedVCF$"#CHROM","_P", imputedVCF$POS,sep="")
   a
}
```

Read in and format the true genotypes
```
true=fread(paste(maindir,"TrueGenotypes.vcf",sep=""),skip=5765)
true2=as.matrix(true[,-c(1:9)])
colnames(true2) = colnames(true)[10:ncol(true)]
rownames(true2) = paste("C",true$"#CHROM","_P",true$POS,sep="")
```


Change "imputed_file1" or "imputed_file2" in the code below to compare accuracies between the different reference populations.

```
imputed_file1="ImputedData_Other_WGS.vcf.gz"
imputed_file2="ImputedData_NZ_WGS.vcf.gz"

imputed1=fread(paste(mkdir,imputed_file1,sep=""))
imputed2=fread(paste(mkdir,imputed_file2,sep=""))

genos1 = getGenotypeCalls(imputed1) 
genos2 = getGenotypeCalls(imputed2) 

dosage1 = getDosages(imputed1) 
dosage2 = getDosages(imputed2) 


##Concordance##
genoConc1 = genos1 == true2  
genoConc12 = genoConc1  
genoConc12[which(is.na(true2))] = NA
IndAcc1 = colMeans(genoConc12,na.rm=T) ## Individual Concordance
summary(IndAcc1)
SNPAcc1 = rowMeans(genoConc12,na.rm=T) ## SNP Concordance
summary(SNPAcc1)

genoConc2 = genos2 == true2  
genoConc22 = genoConc2  
genoConc22[which(is.na(true2))] = NA
IndAcc2 = colMeans(genoConc22,na.rm=T) ## Individual Concordance
summary(IndAcc2)
SNPAcc2 = rowMeans(genoConc22,na.rm=T) ## SNP Concordance
summary(SNPAcc2)


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

#### Extra Exploration of Data
Merge the imputed_files with the info file and you can further explore:
- Does the imputation accuracy depend on the breed? Is the same reference set the best for all breeds? You can use the aggregate function to look at this.
- Try plotting the PCA plot from above but change the size of the points to reflect imputation accuracy. Are there still animals that are poorly imputed?
