# OVERVIEW
In imputation studies, there are two datasets involved: a reference and a target. The reference dataset is a set of individuals that have genotype data at the desired density, while the target data includes individuals that you wish to infer the higher density genotypes. A crucial aspect to any imputation study is identifying the appropriate reference set of individuals to represent the target individuals. However, scenario where only a limited subset of individuals can be genotyped (e.g. due to cost) at the desired density, choosing the right individuals to serve as a reference set become critically important.

# AIMS:
1)	Assess the impact that using different individuals in a reference set has on imputation accuracy.

# BACKGROUND
Due to time and computational limitations, this practical is focused on using data that has already been imputed for you. Briefly, there were x steps: 
1) Genotypes were phased
2) 500kb haplotypes were generated
3) LPChoose was run to identify a set of animals that would capture additional haplotypes within the population (on top of those captured by WGS animals)
4) A variety of reference populations were generated based on country of origin of the WGS animals, LPChoose selections and Randomly selected animals
5) Imputation was run for Chromosome 26 in BEAGLE with the different reference sets
- Effective population size parameter set to 500
- 52 SNPs from the 50k set of SNPs were masked for all ~36k target animals to test imputation accuracy
6) True and imputed genotypes for the 52 SNPs and ~36k target animals were subsetted from the full dataset


# TOOLS:
R: Packages data.table and reshape2 – offer quicker manipulations of large data frames

# INPUT FILES:
- infoFile.csv: Breed information on each animal as well as Principal Components from full genome
- TrueGenotypes.vcf: True Genotypes on each animal to compute accuracy statistics
- ImputedData_xxx.vcf: Imputed genotypes based on reference population xxx with options
   - All WGS animals: ImputedData_ALL_WGS.vcf
   - NZ WGS animals: ImputedData_NZ_WGS.vcf
   - NZ + AUS WGS animals: ImputedData_NZ_AUS_WGS.vcf
   - AUS WGS animals: ImputedData_AUS_WGS.vcf
   - Other WGS animals: ImputedData_Other_WGS.vcf
   - Other + AUS animals: ImputedData_Other_AUS_WGS.vcf
   - LPChoose + WGS animals: ImputedData_LPChoose_WGS.vcf
   - Random + WGS animals: ImputedData_Random_WGS.vcf   
   
   
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
maindir = "/nesi/nobackup/nesi02659/SEP28/practical1/"
```

Read in the info file and take a look at it
```
infoFile = read.csv(paste(maindir,"infoFile.csv",sep=""),header=T)
head(infoFile)
 
```

A PCA plot can give us an idea of how closely related the reference and target animals are. I have already generated the PCs using the full genome (rather than the single chromosome used in this example) and this was done using the "eigen" command run on the genomic relationship matrix. 

The "PlotColour" column has been set for different colours based on breed, with WGS animals a different colour, based on either "NZ" or "Other" as the country that provided the WGS, and LPChoose or Randomly selected animals different colours.
```
par(bg = "gray") # set a gray background so you can see the colours more clearly
plot(infoFile$PC1,infoFile$PC2,pch=19,col=as.character(infoFile$PlotColour),xlab="PC1",ylab="PC2")
points(infoFile$PC1[infoFile$LPChoose],infoFile$PC2[infoFile$LPChoose],pch=19,col=as.character(infoFile$PlotColour[infoFile$LPChoose]),cex=2)
points(infoFile$PC1[infoFile$Random],infoFile$PC2[infoFile$Random],pch=19,col=as.character(infoFile$PlotColour[infoFile$Random]),cex=2)
legend("topright",c("NZ_WGS","Other_WGS","LPChoose","Random"),pch=19,col=c("cyan","blue","green","red"))
```

Now that you have the PCA plot you can look at the relationships between those animals that have WGS and will be used as the reference set of animals, and the target population.

#### Extra Exploration of Data
Some other things you can try by changing the "PlotColour" column in the above code"
- How do Australian WGS animals compare with Other WGS animals?
- How representative are the samples that were selected by LPChoose (LPChoose column)? These were selected based on the WGS animals already being in the reference population.

## 3) Impact of Different Imputation Sets
Here, we’ll test the impact of using different reference sets for imputation of individuals with 50k genotypes to HD genotypes. We have run BEAGLE 5.1 a number of times using a variety of reference populations and provided the output files for you to evaluate imputation accuracy. These files contain the imputed genotypes for the target animals (i.e. animals with 50k genotypes) for a set of SNPs that have been masked to test imputation accuracy.

Reference Population:
- All WGS animals: ImputedData_ALL_WGS.vcf
- NZ WGS animals: ImputedData_NZ_WGS.vcf
- NZ + AUS WGS animals: ImputedData_NZ_AUS_WGS.vcf
- AUS WGS animals: ImputedData_AUS_WGS.vcf
- Other WGS animals: ImputedData_Other_WGS.vcf
- Other + AUS animals: ImputedData_Other_AUS_WGS.vcf
- LPChoose + WGS animals: ImputedData_LPChoose_WGS.vcf
- Random + WGS animals: ImputedData_Random_WGS.vcf

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

getAccuracy = function(test, true,by="row"){
   test=apply(test,2,function(d){as.numeric(as.character(d))})
   true=apply(true,2,function(d){as.numeric(as.character(d))})
   acc = c()
   if(by == "row"){
      for(i in 1:nrow(test)){
         acc = c(acc,cor(test[i,],true[i,],use="complete.obs"))
      }
   }else{
      for(i in 1:ncol(test)){
         acc = c(acc,cor(test[,i],true[,i],use="complete.obs"))
      }
   }
   acc
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
imputed_file1="ImputedData_Other_WGS.vcf"
imputed_file2="ImputedData_NZ_WGS.vcf"

imputed1=fread(paste(maindir,imputed_file1,sep=""))
imputed2=fread(paste(maindir,imputed_file2,sep=""))

genos1 = getGenotypeCalls(imputed1) 
genos1[1:5,1:5]
genos2 = getGenotypeCalls(imputed2) 
genos2[1:5,1:5]

dosage1 = getDosages(imputed1) 
dosage1[1:5,1:5]
dosage2 = getDosages(imputed2) 
dosage2[1:5,1:5]
```

Get the Accuracy of the Imputed Genotype Calls based on concordance with true genotypes.
- IndAcc1 = Accuracy of individuals from imputed_file1
- SNPAcc1 = Accuracy of SNPs from imputed_file1
- IndAcc2 = Accuracy of individuals from imputed_file2
- SNPAcc2 = Accuracy of SNPs from imputed_file2

```
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
```

Get Dosage R2 estimates from BEAGLE output. These are an internal estimate of SNP accuracy. How do these compare with the true accuracy from above?
- dr1 is the dosage R2 for imputed_file1 and can be compared to SNPAcc1
- dr2 is the dosage R2 for imputed_file2 and can be compared to SNPAcc2
```
dr1=sapply(strsplit(as.character(imputed1$INFO),";"),"[",1)
dr1=sapply(strsplit(dr1,"="),"[",2)
dr1=as.numeric(dr1)
summary(dr1)

dr2=sapply(strsplit(as.character(imputed2$INFO),";"),"[",1)
dr2=sapply(strsplit(dr2,"="),"[",2)
dr2=as.numeric(dr2)
summary(dr2)
```

Another way to calculate accuracy is to get the correlation between the true and imputed genotypes or dosages. This method is generally more similar to the R2 estimates from BEAGLE.

SNP accuracy from genotype calls
```
GenoSNPCor1 = getAccuracy(genos1,true2,"row")
summary(GenoSNPCor1)

GenoSNPCor2 = getAccuracy(genos2,true2,"row")
summary(GenoSNPCor2)
```

Individual accuracy from genotype calls
```
GenoIndCor1 = getAccuracy(genos1,true2,"col")
summary(GenoIndCor1)

GenoIndCor2 = getAccuracy(genos2,true2,"col")
summary(GenoIndCor2)
```

SNP accuracy from dosage calls
```
DosSNPCor1 = getAccuracy(dosage1,true2,"row")
summary(DosSNPCor1)

DosSNPCor2 = getAccuracy(dosage2,true2,"row")
summary(DosSNPCor2)
```

Individual accuracy from dosage calls
```
DosIndCor1 = getAccuracy(dosage1,true2,"col")
summary(DosIndCor1)

DosIndCor2 = getAccuracy(dosage2,true2,"col")
summary(DosIndCor2)
```

Which reference set gave the most accurate imputation accuracies?
- DosIndCor1 is the individual accuracy from dosage calls using imputed_file1
- DosIndCor2 is the individual accuracy from dosage calls using imputed_file2

```
t.test(DosIndCor1, DosIndCor2, paired = TRUE, alternative = "two.sided")
```

#### Extra Exploration of Data
Merge the imputed_files with the info file and you can further explore:
- Does the imputation accuracy depend on the breed? Is the same reference set the best for all breeds? You can use the aggregate function to look at this.
- Try plotting the PCA plot from above but change the size of the points to reflect imputation accuracy. Are there still animals that are poorly imputed?
