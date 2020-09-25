OVERVIEW: In imputation studies, there are two datasets involved: a reference and a target. The reference dataset is a set of individuals that have genotype data at the desired density, while the target data includes individuals that you wish to infer the higher density genotypes. A crucial aspect to any imputation study is identifying the appropriate reference set of individuals to represent the target individuals. However, scenario where only a limited subset of individuals can be genotyped (e.g. due to cost) at the desired density, choosing the right individuals to serve as a reference set become critically important.

AIMS:
1)	Give an overview of a haplotype-based method to prioritize selection of reference animals.
2)	Assess the impact that using different individuals in a reference set has on imputation accuracy.

TOOLS:
BCFtools: basic bioinformatics software, in this tutorial, we use it for creating subsets and quality control (http://samtools.github.io/bcftools/bcftools.html)

R: Packages data.table and reshape2 – offer quicker manipulations of large data frames

Beagle 5.0: software for phasing and imputation (https://faculty.washington.edu/browning/beagle/beagle.html)

Julia: Package LPChoose – Code to select individuals based on haplotypes
(https://github.com/reworkhow/LPChoose.jl)

Plinkv1.9: Offers many genetics tools. Simply used to convert vcf to 0/1/2 genotypes here
(https://www.cog-genomics.org/plink2)

cd /nesi/nobackup/nesi02659/SEP28

module load BCFtools
module load VCFtools
module load R
module load Beagle/5.0-12Jul19.0df
module load Julia


FILES:
HD_50K_Overlap_SNPs_26_All_Ans.vcf.gz

STEP 1: Phase all individuals using a common set of genotypes to obtain haplotypes

Background: In this dataset, we have individuals SNP data from whole-genome sequence, a high-density (600k) genotype array, and a medium-density (50k) genotype array. A common set of SNPs were used to combine the individuals into a single set. This common set subsequently underwent QC, and the resulting SNPs were kept for downstream analysis. An example of this is the provided file HD_50K_Overlap_SNPs_26_All_Ans.vcf.gz. 

Action: We will phase this set of individuals using Beagle5.0: 
beagle gt=HD_50K_Overlap_SNPs_26_All_Ans.vcf.gz nthreads=1 ne=500 out=HD_50K_Overlap_SNPs_26_All_Ans_Phased

STEP 2: Convert haplotypes into haploblocks

Background: We can then convert this into haplotypes of a defined size (called a haploblock). Here, we’ll use a haploblock size of 500kb. Note: This is code for a single chromosome – this would be done for all chromosomes and combined into a single file.

Action: This will be achieved in R:
R
require(data.table)
c=26
a=fread(paste("HD_50K_Overlap_SNPs_",c,"_All_Ans_Phased.vcf.gz",sep=""))
win=5e5
b=a[which(a$POS<=win),]

library(reshape2)
x2= b[,-c(1:9)]
x2=as.data.frame(unlist(x2))
colnames(x2)="var"
x2 <- colsplit(x2$var, "[|]", c("left", "right"))
p1=t(matrix(x2$left,nrow(b),ncol(b)-9))
p2=t(matrix(x2$right,nrow(b),ncol(b)-9))

p1a=apply(p1,1,function(x){paste(x,collapse="")})
p2a=apply(p2,1,function(x){paste(x,collapse="")})
p=cbind(p1a,p2a)
pb=as.numeric(as.factor(p))
pb=matrix(pb,nrow(p),ncol(p))

for(i in 2:ceiling(max(a$POS)/win)){
	print(paste(i,"out of",ceiling(max(a$POS)/win),"windows"))
	b=a[which(a$POS>(i-1)*win & a$POS<=i*win),]
	x2= b[,-c(1:9)]
	x2=as.data.frame(unlist(x2))
	colnames(x2)="var"
	x2 <- colsplit(x2$var, "[|]", c("left", "right"))
	p1=t(matrix(x2$left,nrow(b),ncol(b)-9))
	p2=t(matrix(x2$right,nrow(b),ncol(b)-9))

	p1a=apply(p1,1,function(x){paste(x,collapse="")})
	p2a=apply(p2,1,function(x){paste(x,collapse="")})
	p=cbind(p1a,p2a)
	pbs=as.numeric(as.factor(p))
	pbs=matrix(pbs,nrow(p),ncol(p))

	pb=cbind(pb,pbs)
}
pbf=cbind(1:nrow(pb),pb)
pbf=as.data.frame(pbf)
pbf[1:10,1:10]
fwrite(pbf,paste("HD_50K_Overlap_SNPs_",c,"_All_Ans.hap500",sep=""),col.names=F,row.names=F,sep="\t")
q()
n

Background: The above R code provides us with rows=number of animals and columns=2*number of haploblocks plus 1. The first column in this matrix is the individual ID, and the elements of the rest of this matrix are the haploblock alleles each diploid individual possesses. As mentioned earlier, the dataset includes all individuals (WGS, HD, and 50k). For the purposes of this exercise, we will only consider individuals that have HD or WGS genotypes when considering additional animals to include in the reference set. To do this, we can use a weights file to indicate whether LPChoose should consider the individual as a candidate. Additionally, we have selected a random set of 96 individuals for comparison of the performance of LPChoose:
"LPChoose_weights.txt
"Random96.txt" 

STEP 3: Run LPChoose
Background: Using these weights, we can run LPChoose to identify a set of individuals to select in addition to those we already have for an improved reference set. We can also get a feel for haplotype diversity in the set by looking at the number of iterations to capture all haplotypes.
Action: This will be done in Julia. Note: each time LPChoose is run, it generates 3 files (identified_animals.txt, haplotype_coverage.txt, genome_coverage.txt). The file names will need to be changed before next run of LPChoose to prevent them from being overwritten. Also, the example uses a single chromosome to save time. The next step uses outputs for the whole genome.

julia
##LPChoose code with functions to run algorithm
include("LPChoose_weights.jl")

##Weights file to only consider animals with HD genotypes##
weights=readdlm("LPChoose_weights.txt")
weights=weights[1:36122]'

##Specify animals with WGS data – last 935 individuals in our files##
animals_selected=collect(36123:37057)

##Specify a random set of individuals for comparison##
ran=round.(Int,readdlm("Random96.txt"))
animals_selected_ran=[animals_selected;ran]

##Assess haplotype diversity in whole dataset##
LPChoose("HD_50K_Overlap_SNPs_26_All_Ans.hap500",2000,0.0,nsteps=2000)

##Choose additional 96 individuals using LPChoose algorithm##
LPChoose("HD_50K_Overlap_SNPs_26_All_Ans.hap500",96,0.0,nsteps=96,preselected_animals=animals_selected,weights=weights)

##Assess performance of randomly selected 96 individuals##
LPChoose("HD_50K_Overlap_SNPs_26_All_Ans.hap500",0,0.0,nsteps=1,preselected_animals=animals_selected_ran)
STEP 4: Assess haplotype diversity

Background: We can see how many iterations are required to capture all haplotypes, how representative of our target population our reference is, etc. We’ll use the outputs for all haplotypes.

Action: This can be in R:

##All individuals  – no pre-selection##
x1=read.csv("identified_animals_hap500_all.txt",header=T)
x2=read.csv("haplotype_coverage_hap500_all.txt",header=T)
x3=read.csv("genome_coverage_hap500_all.txt",header=T)
xx=merge(x1,x2,by="step")
x=merge(xx,x3,by="step")
colnames(x)=c("Step","AnID","PropHap","PropGen")

##Considering individuals with WGS and selecting an additional 96 animals##
u1=read.csv("identified_animals_hap500_wgs.txt",header=T)
u2=read.csv("haplotype_coverage_hap500_wgs.txt",header=T)
u3=read.csv("genome_coverage_hap500_wgs.txt",header=T)
uu=merge(u1,u2,by="step")
u=merge(uu,u3,by="step")
colnames(u)=c("Step","AnID","PropHap","PropGen")

##Considering animals with WGS and additional randomly selected 96 individuals##
r1=read.csv("identified_animals_hap500_ran.txt",header=T)
r2=read.csv("haplotype_coverage_hap500_ran.txt",header=T)
r3=read.csv("genome_coverage_hap500_ran.txt",header=T)
rr=merge(r1,r2,by="step")
r=merge(rr,r3,by="step")
colnames(r)=c("Step","AnID","PropHap","PropGen")

jpeg("PropHaps_Full_500.jpg",height=800,width=800)
botleft <- c(0.000001,-0.0400000)
topright <- c(10^4.5339724,1.0400000)
par(cex.lab=1.5)
plot(c(1:22887), c(1:22887)/22887,log="x", type = "n", axes = FALSE, xlab="Step", ylab="Proportion")
lim <- par("usr")
rect(botleft[1], botleft[2], topright[1], topright[2], border = "grey",
     col = "grey")
axis(1,cex.axis=1.5) ## add axes back
axis(2,cex.axis=1.5)
box()  
points(x$PropHap~x$Step,col=1,type="l",lwd=2)
points(x$PropGen~x$Step,col=1,type="l",xlim=c(1,22887),lwd=2,lty=2)
segments(x0=100,y0=0,x1=100,y1=x[101,"PropGen"],col=1,lty=2)
segments(x0=1,y0=x[101,"PropGen"],x1=100,y1=x[101,"PropGen"],col=1,lty=2)
segments(x0=500,y0=0,x1=500,y1=x[501,"PropGen"],col=2,lty=2)
segments(x0=1,y0=x[501,"PropGen"],x1=501,y1=x[501,"PropGen"],col=2,lty=2)
segments(x0=1000,y0=0,x1=1000,y1=x[1001,"PropGen"],col=4,lty=2)
segments(x0=1,y0=x[1001,"PropGen"],x1=1000,y1=x[1001,"PropGen"],col=4,lty=2)
segments(x0=100,y0=0,x1=100,y1=x[101,"PropHap"],col=1)
segments(x0=1,y0=x[101,"PropHap"],x1=101,y1=x[101,"PropHap"],col=1)
segments(x0=500,y0=0,x1=500,y1=x[501,"PropHap"],col=2)
segments(x0=1,y0=x[501,"PropHap"],x1=501,y1=x[501,"PropHap"],col=2)
segments(x0=1000,y0=0,x1=1000,y1=x[1001,"PropHap"],col=4)
segments(x0=1,y0=x[1001,"PropHap"],x1=1000,y1=x[1001,"PropHap"],col=4)
dev.off()

jpeg("PropHaps_Full_500_WGS.jpg",height=800,width=800)
botleft <- c(0.000001,-0.0400000)
topright <- c(10^4.5339724,1.0400000)
par(cex.lab=1.5)
plot(c(1:22887), c(1:22887)/22887,log="x", type = "n", axes = FALSE, xlab="Step", ylab="Proportion")
lim <- par("usr")
rect(botleft[1], botleft[2], topright[1], topright[2], border = "grey",
     col = "grey")
axis(1,cex.axis=1.5) ## add axes back
axis(2,cex.axis=1.5)
box()  
points(x$PropHap~x$Step,col=1,type="l",lwd=2)
points(x$PropGen~x$Step,col=1,type="l",lwd=2,lty=2)
abline(h=u[1,"PropGen"],lty=2,col="blue",lwd=2)
abline(h=u[1,"PropHap"],col="blue",lwd=2)
abline(h=u[nrow(u),"PropGen"],lty=2,col="cyan",lwd=2)
abline(h=u[nrow(u),"PropHap"],col="cyan",lwd=2)
abline(h=r[nrow(r),"PropGen"],lty=2,col="red",lwd=2)
abline(h=r[nrow(r),"PropHap"],col="red",lwd=2)
dev.off()

STEP 5: Look at population structure

Background: We can also have a look at the population structure through a PCA, and see where on the individuals selected are in the PCA. Note: I have already generated the PCs, but example code is provided for how to do this.

##LOOK AT POPULATION STRUCTURE THROUGH PCA PLOT - SINGLE CHROMOSOME EXAMPLE###
##CONVERT VCF TO 012 FORMAT FOR PCA##
##GRM AND PRINCIPAL COMPONENTS##
#conda activate plink
#plink --vcf HD_50K_Overlap_SNPs_26_All_Ans.vcf.gz --recodeA --out HD_50K_Overlap_SNPs_26_All_Ans.012
#R
#require(data.table)
#c=26
#x=fread(paste("HD_50K_Overlap_SNPs_",c,"_All_Ans.012.raw",sep=""))
#require(zoo)
#x=x[,-c(1:6)]
#x=na.aggregate(x)
#p=colMeans(x)/2
#x=scale(x, scale = FALSE)
#x=(tcrossprod(x))
#G=x/sum(2*p*(1-p))
#pca=eigen(G)
#pc=pca$vectors
#pr=pca$values/sum(pca$values)
##

Action: This can be in R:
R
require(data.table)

##Read in Principal Components##
pc=fread("10PCs_HD_50K_Overlap_with_metadata.pc") 
pr=fread("10PCs_HD_50K_Overlap.pr") 

##INCLUDE ADDITIONAL CANDIDATES##
lpc=fread("LPChoose_HD.ans",header=F)
lpc=lpc[936:nrow(lpc),]
rand=fread("Random_HD.ans",header=F)
rand=rand[1:96,]

jpeg("PCA.jpg",height=1200,width=2400)
par(cex.lab=1.5)
plot(x=c(min(pc[,PC1]),max(pc[,PC1])),y=c(min(pc[,PC2]),max(pc[,PC2])), type = "n", axes = FALSE, xlab=paste("PC1 (",round(pr[1,]*100,1),"%)",sep=""), ylab=paste("PC2 (",round(pr[2,]*100,1),"%)",sep=""))
lim <- par("usr")
rect(-0.01274953, -0.01701290, 0.02107473, 0.03045960, border = "grey",col = "grey")
axis(1,cex.axis=1.5) ## add axes back
axis(2,cex.axis=1.5)
box()  
points(pc[,PC2]~pc[,PC1],pch=19)
points(pc[which(pc$ChipBreed=="Texel"),PC2]~pc[which(pc$ChipBreed =="Texel"), PC1],col=7,pch=19)

points(pc[which(pc$ChipBreed =="FinnX"), PC2]~pc[which(pc$ChipBreed =="FinnX"), PC1],col="orange",pch=19)

points(pc[which(pc$ChipBreed =="Primera"), PC2]~pc[which(pc$ChipBreed =="Primera"), PC1],col="pink",pch=19)

points(pc[which(pc$ChipBreed =="Peren"), PC2]~pc[which(pc$ChipBreed =="Peren"), PC1],col=6,pch=19)

points(pc[which(pc$ChipBreed =="Rom"), PC2]~pc[which(pc$ChipBreed =="Rom"), PC1],col="firebrick4",pch=19)

points(pc[which(pc$ChipBreed =="Coop"), PC2]~pc[which(pc$ChipBreed =="Coop"), PC1],col="white",pch=19)

points(pc[which(pc$WGSCountry!=""), PC2]~pc[which(pc$WGSCountry!=""), PC1],col="blue",pch=19,cex=2)

points(pc[which(pc$ WGSCountry =="NZ"),PC2]~pc[which(pc$ WGSCountry =="NZ"),PC1],col="cyan",pch=19,cex=2)

points(pc[which(ind2$V1 %in% lpc$V1),pcy]~pc[which(ind2$V1 %in% lpc$V1),pcx],col="green",pch=19,cex=2)

points(pc[which(ind2$V1 %in% rand$V1),pcy]~pc[which(ind2$V1 %in% rand$V1),pcx],col="red",pch=19,cex=2)

dev.off()
q()
n

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

##ASSESS IMPUTATION ACCURACY FOR THE DIFFERENT SETS - CHANGE "imputed_file"###
R
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
