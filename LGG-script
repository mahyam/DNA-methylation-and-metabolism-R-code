# READING DNA METHYLATION FILE NAMES
 
setwd("/Volumes/X/LGG 2 2/DNA_Methylation/1")
filenamesmeth<- list.files("LGG METHYLATION")
# READING GENE EXPRESSION FILE NAMES

setwd("/Volumes/X/LGG 2 2/RNASeqV2/1")
filenamesexp<- list.files("LGG EXPRESSION")

#USING TCGA BARCODES TO FIND AND MATCH SAMPLES THAT HAVE BOTH METHYLATION AND EXPRESSION DATA

setwd("/Volumes/X/LGG 2 2")
MAP<- read.table("FILE_SAMPLE_MAP.txt",header=TRUE)

rnaseqind=c()
for(i in 1:length(filenamesexp)){
	rnaseqind[i]<- which(MAP$filename==filenamesexp[i])
}

rnaseqbars<- MAP$barcode.s.[rnaseqind]
expbars<- rnaseqbars

# BARCODES:

m<- length(filenamesmeth)
e<- length(expbars)
s<- max(m,e)

barcodes<- matrix(0,s,2)

for(i in 1:m){
for(j in 1:e){
if(substr(expbars[j],1,15) == substr(filenamesmeth[i],46,60)){
barcodes[i,1]<- i
barcodes[i,2]<- j
}
else if(substr(expbars[j],1,15) == substr(filenamesmeth[i],45,59)){
barcodes[i,1]<- i
barcodes[i,2]<- j
}

}}

NMETH<- filenamesmeth[barcodes[,1]]
NEXP<- expbars[barcodes[,2]]

setwd("/Volumes/X/LGG 2 2/")
write.table(NEXP, "LGG BARCODES 534.txt")


# READING GENE EXPRESSION DATA 

setwd("/Volumes/X/LGG 2 2/RNASeqV2/1/LGG EXPRESSION")

expind<- filenamesexp[barcodes[,2]]
N<- length(expind)
LGGEXP=matrix(0,20531,N)
for(i in 1:N){
N2<- read.table(expind[i],header=TRUE)
LGGEXP[,i]<- as.numeric(N2[,2])
}

setwd("/Volumes/X/LGG 2 2")

write.table(LGGEXP,"LGG expression of 534.txt")
write.table(N2[,1],"LGG expression gene list.txt")


 #**********#
 # FIGURE 1 #
 #**********#
 
# READING DNA METHYLATION DATA  

setwd("/home/X/LGG")
filenamesmeth<- list.files("LGG METHYLATION")
barcodes<- read.table("LGG BARCODES 534.txt")
N=dim(barcodes)[1]

inds<- matrix(0,N,2)
for(i in 1:N){
	inds[i,1]<- i
for(j in 1:length(filenamesmeth)){
if(substr(barcodes[i,1], 1,15) == substr(filenamesmeth[j],46,60)){
inds[i,2]<-j
}
else if(substr(barcodes[i,1],1,15) == substr(filenamesmeth[j],45,59)){
inds[i,2]<-j
}
else if(substr(barcodes[i,1],1,15) == substr(filenamesmeth[j],47,61)){
inds[i,2]<-j
}

}}
inds


setwd("/home/X/LGG/LGG METHYLATION")

METHYLATION=matrix(0,485577,(dim(inds)[1]))
for(i in 1:(dim(inds)[1])){
T<- read.table(filenamesmeth[inds[i,2]],header=F,stringsAsFactors=F,skip=2,sep="\t")
METHYLATION[,i]<- T[,2]
cat(i)
}


# Filter out probes with more than 80% missing values:

N<- dim(METHYLATION)[2]
badp=c()
for(i in 1:485577){
	if(length(which(is.na(METHYLATION[i,])==TRUE))>(N*0.8)  ){
		badp<- c(badp,i)
	}
}

METHYLATION<- METHYLATION[-badp,]

# measure global methylation:

averagemeth=c()
for(i in 1:N){
	averagemeth[i]<- mean(METHYLATION[,i], na.rm=TRUE)
}

range(averagemeth)


setwd("/home/X/LGG")
write.table(averagemeth,"LGG average global meth after filtering bad probes.txt")


# reading Illumina 450k chip probe annotation information

table<- read.table("meth probe info2.txt",header=TRUE)
table<- table[-badp,]
allanot<- T[-badp,4:5]
colnames(allanot)<- c("chr","position")

#  EXCLUDE SEX CHROMOSOMES AND SORT THE METHYLATION MATRIX:

X<- which(table$Chromosome=="X")
Y<- which(table$Chromosome=="Y")
sex<- c(X,Y)

TAB<- data.frame(allanot[-sex,],table[-sex,],METHYLATION[-sex,])

TAB<- TAB[with(TAB,order(chr,position)),]

# removing extra chromosomes

unique(TAB$chr)

NN<- which(is.na(TAB$chr)==TRUE)

TAB<- TAB[-NN,]
XX<- which(TAB$chr=="X")

TAB<- TAB[-XX,]

unique(TAB$chr)

# LOCAL DNA METHYLATIONS =10 KB BINS
 # average window analysis

chroms<- as.numeric(unique(TAB$chr))
bins<- matrix(0,1,(N+2))
colnames(bins)<- c("chr","range",seq(1:N))
for(i in 1:length(chroms)) {
	ind<- which(TAB$chr==chroms[i])
	tab<- TAB[ind,]
start<-tab$position[1]
end<- tab$position[dim(tab)[1]]
b<- seq(start,end,by=10000)
bins<- rbind(bins,c(chroms[i],0,rep(0,N)))
s<- seq(start,end,by=10000)
s<- c(s,end)

ave<-matrix(0,length(b),N)
        	for(j in 1:N){
    		ave[,j]<- tapply(tab[,j+11], cut(tab$position,s , include.lowest=TRUE), mean,na.rm=TRUE)

    }
range<- levels(cut(tab$position, s, include.lowest=TRUE))    
ch<- rep(chroms[i],length(b))	 	 	
newbin<- cbind(ch,range,ave)
bins<- rbind(bins,newbin)
	 			 	
	 cat(i)
}
setwd("/home/X/LGG")
write.table(bins, "all 10kb bin ave meths in LGG.txt")

dim(bins)

# filter out non-variable bins (sd<0.2)

bins[which(bins=="NaN")]<- NA
bad<- c()
for( i in 1:dim(bins)[1]){
	if(length(which(is.na(bins[i,])==TRUE))>(N*0.9)){
		bad<- c(bad,i)
	}
	
	else if(sd(as.numeric(bins[i,3:(N+2)]),na.rm=TRUE)<0.2){
		bad<- c(bad,i)
	}
}

length(bad)

response<- bins[-bad,]
dim(response)

write.table(response, "LOCAL methylation matrix in LGG.txt")

r<- matrix(as.numeric(unlist(response[,3:(N+2)])),nrow=nrow(response))

box<- colSums(r, na.rm=TRUE)/(dim(r)[1])
range(box)

write.table(box,"average local meth for boxplots in lgg.txt")


# Correlation with methionine cycle genes

setwd("/home/X/LGG")
LGGEXP<- data.matrix(read.table("LGG expression of 534.txt"))

# filter out low expression genes based on RSEM 3

LGGEXP[which(LGGEXP==0)]<- 1
lLGGEXP<- log2(LGGEXP)
N=dim(LGGEXP)[2]
nbads=c()

for(i in 1: (dim(lLGGEXP)[1])){
	if(length(which(lLGGEXP[i,]>=3))<(N*0.7)){
		nbads<- c(nbads,i)	
	}
}

nLGGEXP<- lLGGEXP[-nbads,]

# read gene names as ordered in TCGA's RNA-seq data

genenames<- read.table("expression gene list.txt", header=TRUE)

nGenenames<- genenames[-nbads,1]

MECEXP<- c(which(nGenenames=="MAT2B"), which(nGenenames=="MTR"), which(nGenenames=="BHMT2"), which(nGenenames=="AHCY"))
nGenenames[MECEXP]
# MAT2B MTR   BHMT2 AHCY

  # correlated Met cycle gene expression with global methylation:
  
globalC=c()
globalP=c()
for(i in 1:(dim(nLGGEXP)[1])){
	globalC[i]<- cor.test(nLGGEXP[i,],averagemeth,method="spearman")$estimate
	globalP[i]<- cor.test(nLGGEXP[i,],averagemeth,method="spearman")$p.value

	
}

range(globalC)


-log10(globalP[MECEXP])

globalC[MECEXP]


setwd("/home/X/LGG")
write.table(cbind(globalP,globalC), "genome wide LGG cors with average global meth.txt")

  #correlated met cycle gene expression with local DNA methylation  (median pval)

  
MAT2B<- nLGGEXP[which(nGenenames=="MAT2B"),]
MAT2BCORS=c()
	for (i in 1:(dim(response)[1])){
		MAT2BCORS[i]<- cor.test(MAT2B,as.numeric(response[i,3:(N+2)]), method="spearman" )$p.value

	}
median(MAT2BCORS,na.rm=TRUE)

 MAT2Bl<-  (length(which(MAT2BCORS<0.01))*100)/(dim(response)[1])
  MAT2Bl

MTR<- nLGGEXP[which(nGenenames=="MTR"),]
MTRCORS=c()
	for (i in 1:(dim(response)[1])){
		MTRCORS[i]<- cor.test(MTR,as.numeric(response[i,3:(N+2)]), method="spearman" )$p.value

	}
median(MTRCORS,na.rm=TRUE)
MTRl<-  (length(which(MTRCORS<0.01))*100)/(dim(response)[1])
MTRl



AHCY<- nLGGEXP[which(nGenenames=="AHCY"),]
AHCYCORS=c()
	for (i in 1:(dim(response)[1])){
		AHCYCORS[i]<- cor.test(AHCY,as.numeric(response[i,3:(N+2)]), method="spearman" )$p.value

	}

median(AHCYCORS,na.rm=TRUE)

AHl<- (length(which(AHCYCORS<0.01))*100)/(dim(response)[1])
AHl



BHMT2<- nLGGEXP[which(nGenenames=="BHMT2"),]
BHCORS=c()
	for (i in 1:(dim(response)[1])){
		BHCORS[i]<- cor.test(BHMT2,as.numeric(response[i,3:(N+2)]), method="spearman" )$p.value

	}

median(BHCORS,na.rm=TRUE)


BHl<- (length(which(BHCORS<0.01))*100)/(dim(response)[1])
BHl

# Supp F1 #

# Global cors vs. local cors for 500 RANDOM genes in the genome :

random<- sample(seq(1,(dim(nLGGEXP)[1])),500)
randgenes<- nGenenames[random]
   
RANDCORS=matrix(0,500,(dim(response)[1]))  

for(i in 1:500){
	for (j in 1:(dim(response)[1])){
		RANDCORS[i,j]<- cor.test(nLGGEXP[random[i],],as.numeric(response[j,3:(N+2)]), method="spearman" )$p.value
	}
	cat(i)
} 
       
    
sigcors=c()
for(i in 1:500){
	sigcors[i]<- (length(which(RANDCORS[i,]<0.01))*100)/(dim(response)[1])
    
} 

A<- mean(sigcors,na.rm=TRUE)
A


GRANDCORS=c()  

	for (i in 1:500){
		GRANDCORS[i]<- cor.test(nLGGEXP[random[i],],averagemeth, method="spearman" )$p.value
	}
	
B<- (length(which(GRANDCORS<0.01))*100)/500
B

fisher.test(matrix(c(length(which(RANDCORS<0.01)),(500*297),length(which(GRANDCORS<0.01)),500),ncol=2))


range(sigcors)

write.table(sigcors, "random 500 genes percent of sig cors with local meth in lgg.txt")

length(which(sigcors>MAT2Bl))/500

length(which(sigcors>MTRl))/500

length(which(sigcors>AHl))/500

length(which(sigcors>BHl))/500


#**********#
# FIGURE 2 #
# *********#


dim(response)

ranges<- response[,1:2]
response<- response[,-c(1,2)]

df.feat<- data.frame(t(nLGGEXP[MECEXP,]))
colnames(df.feat)<- as.character(nGenenames[MECEXP])

# ***** MODELING ****** #

library("randomForest")
library("glmnet")


# ONLY USING MEC VARIABLES

 # setting number of folds for cv

k=7

test.ind.mat <- matrix(sample(1:dim(response)[2]),nrow=k)

rfMSE.mat <- matrix(0,nrow=dim(response)[1],ncol=k)
elMSE.mat <- matrix(0,nrow=dim(response)[1],ncol=k)
varimp<- matrix(0,1,(dim(df.feat)[2]))
varcoef<- matrix(0,1,((dim(df.feat)[2])+1))

for (i in 1:k) {
	
	# Run random forest and elastic net on each response
	for (j in 1:dim(response)[1]) {
		
		df.resp <- response[j,]
df.full <- data.frame("resp"=as.numeric(df.resp),df.feat)

  # Create train and test set
	df.train <- df.full[-test.ind.mat[i,],]
	df.test <- df.full[test.ind.mat[i,],]
	
	# use na.rough to impute missing values
	df.train.rough <- na.roughfix(df.train)
	df.test.rough <- na.roughfix(df.test)
	

   #rf
		rf <- randomForest(y=df.train.rough[,1],x=df.train.rough[,-1],ntree=500,mtry=3,importance=TRUE)
		rf.pred <- predict(rf,df.test.rough)
		rfMSE.mat[j,i] <- (1/(dim(test.ind.mat)[2]))*(sum((rf.pred - df.test.rough[,1])^2))
		varimp<- rbind(varimp,importance(rf,scale=TRUE)[,1])

    #el
    

Y<- df.train.rough[,1]

Ecvfit = cv.glmnet(data.matrix(df.train.rough[,-1]),Y,alpha=0.5,nfold=5)


tfit<- predict(Ecvfit,newx=data.matrix(df.test.rough[,-1]) ,s="lambda.min")

elMSE.mat[j,i] <- (1/(dim(test.ind.mat)[2]))*(sum((tfit - df.test.rough[,1])^2))
varcoef<- rbind(varcoef,coef(Ecvfit, s = "lambda.min")[,1])
				
		cat("Fold:  ",i,"  of ",k,";",               "response ", j ," . \n")
	}
	
}

setwd("/home/X/LGG")
write.table(varimp,"MEC only var imps.txt")
write.table(varcoef,"MEC only var coefs.txt")
write.table(rfMSE.mat,"MEC only RF MSE mat.txt")
write.table(elMSE.mat,"MEC only EL MSE mat.txt")


# errors 
mean(data.matrix(rfMSE.mat))

mean(data.matrix(elMSE.mat))


rfE<- sqrt(rfMSE.mat)
elE<- sqrt(elMSE.mat)

rfError<- rowSums(rfE)/k
elError<- rowSums(elE)/k

# find regions that were predicted with small error 

length(which(rfError<0.2))
length(which(elError<0.2))


goodrf<- which(rfError<0.25)
goodel<- which(elError<0.25)
goodcommon<- intersect(goodrf,goodel)
response<- read.table("LOCAL smaller methylation matrix in LGG.txt")
ranges<- response[,1:2]
common<- ranges[goodcommon,]
dim(common)
setwd("/home/X/LGG")
write.table(common,"common good bins found by rf and el for MEC only predictions in lgg.txt")

# bin annotation 

library(GenomicRanges)

# create peaks from ranges

commonP<- read.table("common good bins found by rf and el for MEC only predictions in lgg.txt")

temp<- commonP$range

temp<- gsub("]", "", temp)
temp<- substring(temp, 2, nchar(temp))

temp2<- matrix(as.numeric(unlist(strsplit(temp,","))),ncol=2,byrow=TRUE)

dim(temp2)

peaks<- data.frame(chr=commonP$chr,start=temp2[,1],end=temp2[,2]) 

peaks$chr<- paste0("chr",peaks$chr)
peaks.gr= GRanges (seqnames=peaks$chr,ranges=IRanges(start=peaks$start,end=peaks$end))

write.table(as.data.frame(peaks.gr),"common good bins IRANGES from MEC only in lgg.txt")

 # ( use GREAT for dist to TSS)


# MEC variable contributions:

vim<- varimp[-1,]
meanImp<- colSums (vim)/(k*(dim(response)[1]))

meanImp
 
 vco<- varcoef[-1,]
meanCo<- colSums (vco)/((dim(response)[1])*k)

meanCo
nco<- data.matrix(vco)
nco[which(nco!=0)]<- 1
varuse<- colSums(nco)

varuse/(k*(dim(response)[1]))


# MODELING LOCAL DNA METHYLATION USING THE INTEGRATIVE APPROACH (Using ALL variables)

# GENERATING DIFFERENT CLASSES OF VARIABLES FROM GENE EXPRESSION

  #SAM metabolizing/epigenetic:

SAMs<- read.table("SAM metabolizing and HM.txt")
SAM<- unique(SAMs[,1])
	
SM=c()
for(i in 1:75){
	if(is.na(match(as.character(SAM[i]), nGenenames))=="FALSE"){
    	SM[length(SM)+1]<- which(is.na(match(nGenenames,as.character(SAM[i])))=="FALSE")
    	
    	}}
SAMvars<- data.frame(t(nLGGEXP[SM,]))
colnames(SAMvars)<- nGenenames[SM]    

setwd("/home/X/LGG")	
write.table(SAMvars,"SAM metabolizing variables in lgg.txt")   	 	
 
 
  # chromatin remodellers:
setwd("/home/X/PRAD")    	
RMs<- read.table("new remodelers and dna methylation machinery.txt")  	
RM=c()
for(i in 1:75){
	if(is.na(match(as.character(RMs[i,1]), nGenenames))=="FALSE"){
    	RM[length(RM)+1]<- which(is.na(match(nGenenames,as.character(RMs[i,1])))=="FALSE")
    	
    	}}
RMvars<- data.frame(t(nLGGEXP[RM,]))
colnames(RMvars)<- nGenenames[RM]    

setwd("/home/X/LGG")	
write.table(RMvars,"new chromatin remodeller and dna meths variables in lgg.txt")  

  	 	
   # SGOC:
setwd("/home/X/PRAD")	
sgoc<- read.table("newnetALD.txt")
sg=c()
for(i in 1:65){
	if(is.na(match(as.character(sgoc[i,1]), nGenenames))=="FALSE"){
    	sg[length(sg)+1]<- which(is.na(match(nGenenames,as.character(sgoc[i,1])))=="FALSE")
    	
    	}}
SGOCvars<- data.frame(t(nLGGEXP[sg,]))
colnames(SGOCvars)<- nGenenames[sg]

setwd("/home/X/LGG")	
write.table(SGOCvars,"SGOC variables in lgg.txt")  


   # TFs:
TFs<- read.table("TF list for lgg.txt")  	
TF=c()
for(i in 1:dim(TFs)[1]){
	if(is.na(match(as.character(TFs[i,1]), nGenenames))=="FALSE"){
    	TF[length(TF)+1]<- which(is.na(match(nGenenames,as.character(TFs[i,1])))=="FALSE")
    	
    	}}
TFvars<- data.frame(t(nLGGEXP[TF,]))
colnames(TFvars)<- nGenenames[TF]    	
write.table(TFvars,"small list transcription factor variables in lgg.txt")   	
 

#LGG TCGA : cnv and mut from cBioPortal and paper

  # mutation freq higher than %5 
  
mut<- read.table("muts list for lgg.txt", header=TRUE)
bars<- read.table("LGG BARCODES 534.txt")
intersect(substr(bars[,1],1,15), mut[-1,1])
muts<- mut[-1,-1]
colnames(muts)<- as.character(unlist(mut[1,-1]))
colnames(muts) <- paste(colnames(muts),"mut", sep = "_")

  # AMP / deletions with freq higher than %15  

cnv<- read.table("cnvs list for lgg.txt", header=TRUE)
intersect(substr(bars[,1],1,15), cnv[-1,1])
cnvs<- cnv[-1,-1]
colnames(cnvs)<- as.character(unlist(cnv[1,-1]))
colnames(cnvs) <- paste(colnames(cnvs),"cn", sep = "_")

barcodes<- cnv[-1,1]

MutMat=matrix(0,N,dim(muts)[2])
colnames(MutMat)<- colnames(muts)
cnvMat=matrix(0,N,dim(cnvs)[2])
colnames(cnvMat)<- colnames(cnvs)
 for(i in 1:N){
 	ind<- which(barcodes==substr(bars[i,1],1,15))
 	cnvMat[i,]<- as.numeric(as.character(unlist(cnvs[ind,])))
 	p<- which(muts[ind,]!="NaN")
 	if(length(p)>0){
 	MutMat[i,p]<- 1 	
 	}
 		
 } 

write.table(MutMat,"final mutations matrix in lgg.txt")
write.table(cnvMat,"final cnvs matrix in lgg.txt")


#****#

setwd("/home/X/LGG")

ranges<- response[,1:2]
response<- data.matrix(response[,-c(1,2)])

# reading all variables for the comprehensive modeling

Clin<- read.csv("small clinical variables for lgg.csv")
Clin<- Clin[,-1]
SAMvars<- read.table("SAM metabolizing variables in lgg.txt")
RMvars<- read.table("new chromatin remodeller and dna meths variables in lgg.txt")
SGOCvars<- read.table("SGOC variables in lgg.txt")
MEC<-SGOCvars[,c(4,6,8,9)]
SGOCvars<- SGOCvars[,-c(4,6,8,9)]

TFvars<- read.table("small list transcription factor variables in lgg.txt")
MutMat<- read.table("final mutations matrix in lgg.txt")
cnvMat<- read.table("final cnvs matrix in lgg.txt")


X<- cbind(Clin,MutMat,cnvMat,MEC, TFvars,RMvars,SAMvars,SGOCvars)


X1<- X[, unique(colnames(X))]
str(X1)
W<- dim(X1)[2]

k=3

test.ind.mat <- matrix(sample(1:dim(response)[2]),nrow=k)

rfMSE.mat <- matrix(0,nrow=dim(response)[1],ncol=k)
elMSE.mat <- matrix(0,nrow=dim(response)[1],ncol=k)
varimp<- matrix(0,1,W)
varcoef<- matrix(0,1,(W+1))

for (i in 1:k) {
	
	# Run random forest and elastic net on each response
	
	for (j in 1:dim(response)[1]) {
		
		df.resp <- response[j,]
df.full <- data.frame("resp"=as.numeric(df.resp),X1)

  # Create train and test set
	df.train <- df.full[-test.ind.mat[i,],]
	df.test <- df.full[test.ind.mat[i,],]
	
	# use na.rough to impute missing values
	df.train.rough <- na.roughfix(df.train)
	df.test.rough <- na.roughfix(df.test)
	

   #rf
		rf <- randomForest(y=df.train.rough[,1],x=df.train.rough[,-1],ntree=500,importance=TRUE)
		rf.pred <- predict(rf,df.test.rough)
		rfMSE.mat[j,i] <- (1/(dim(test.ind.mat)[2]))*(sum((rf.pred - df.test.rough[,1])^2))
		varimp<- rbind(varimp,importance(rf,scale=TRUE)[,1])

    #el
    

Y<- df.train.rough[,1]

Ecvfit = cv.glmnet(data.matrix(df.train.rough[,-1]),Y,alpha=0.5,nfold=5)


tfit<- predict(Ecvfit,newx=data.matrix(df.test.rough[,-1]) ,s="lambda.min")

elMSE.mat[j,i] <- (1/(dim(test.ind.mat)[2]))*(sum((tfit - df.test.rough[,1])^2))
varcoef<- rbind(varcoef,coef(Ecvfit, s = "lambda.min")[,1])
				
		cat("Fold:  ",i,"  of 3;               response ",j," . \n")
	}
	
}

write.table(varimp,"ALL VARS var imps in lgg.txt")
write.table(varcoef,"ALL VARS var coefs in lgg.txt")
write.table(rfMSE.mat,"ALL VARS RF MSE mat in lgg.txt")
write.table(elMSE.mat,"ALL VARS EL MSE mat in lgg.txt")


mean(data.matrix(elMSE.mat),na.rm=TRUE)

mean(data.matrix(rfMSE.mat),na.rm=TRUE)


rfE<- sqrt(rfMSE.mat)
elE<- sqrt(elMSE.mat)

rfError<- rowSums(rfE)/k
elError<- rowSums(elE)/k

# find regions where methylation was predicted with small error

length(which(rfError<0.2))

length(which(elError<0.2))

goodrf<- which(rfError<0.22)
goodel<- which(elError<0.22)
goodcommon<- intersect(goodrf,goodel)

common<- ranges[goodcommon,]
dim(common)

write.table(table(common$chr),"chromosome distributions for common good bins found by rf and el for all vars  predictions in lgg.txt")

write.table(common,"common good bins found by rf and el for all vars  predictions in lgg.txt")

vim<- varimp[-1,]
meanImp<- colSums (vim)/((dim(response)[1])*3)
 

gc<- c(goodcommon,goodcommon+297,goodcommon+(2*297)) 

GmeanImp<- colSums (vim[gc,])/(dim(common)[1]*3)
sort(GmeanImp)
write.table(GmeanImp,"mean variable importances all variables all low mse models in lgg.txt") 
    
vco<- varcoef[-1,]
nco<- data.matrix(vco)

nco[which(nco!=0)]<- 1
varuse<- colSums(nco)

sort(varuse)

write.table(varuse,"variable usage all variables all models in lgg.txt")

Gvaruse<- colSums(nco[gc,])

sort(Gvaruse)

write.table(Gvaruse,"variable usage all variables low mse models in lgg.txt")

# annotation of bins

library(GenomicRanges)

 # create peaks from ranges

temp<- ranges$range[goodcommon]

temp<- gsub("]", "", temp)
temp<- substring(temp, 2, nchar(temp))

temp2<- matrix(as.numeric(unlist(strsplit(temp,","))),ncol=2,byrow=TRUE)

dim(temp2)

peaks<- data.frame(chr=ranges$chr[goodcommon],start=temp2[,1],end=temp2[,2]) 

peaks$chr<- paste0("chr",peaks$chr)
peaks.gr= GRanges (seqnames=peaks$chr,ranges=IRanges(start=peaks$start,end=peaks$end))

write.table(as.data.frame(peaks.gr)," good bins with IRanges for GREAT all vars model in lgg.txt")

#get distance to nearest TSS

 # ( use GREAT for dist to TSS)

# MEC vars contributions:

MeC<- c(29,30,31,32)
GmeanImp[MeC]

    rank(GmeanImp)[MeC]/W

     


Gvaruse[(MeC+1)]

Gvaruse[(MeC+1)]/(dim(peaks)[1]*3)

      
Gvar<- Gvaruse[-1]

rank(Gvar)[MeC]/W 

  AHCY     MAT2B       MTR     BHMT2 
0.2700893 0.8772321 0.2968750 0.4575893

 
# Contribution of variable CLASSES (average)

which(colnames(SAMvars)=="MAT2B") 
SA<- SAMvars[,-40] 

rfimps<- vim[gc,]

CNVcont<- which(is.na(match(colnames(rfimps),colnames(cnvMat)))==FALSE)
Somacont<-which(is.na(match(colnames(rfimps),colnames(MutMat)))==FALSE)
TFcont<-which(is.na(match(colnames(rfimps),colnames(TFvars)))==FALSE)
Mcont<-which(is.na(match(colnames(rfimps),colnames(MEC)))==FALSE)
SGcont<-which(is.na(match(colnames(rfimps),colnames(SGOCvars)))==FALSE)
SAcont<-which(is.na(match(colnames(rfimps),colnames(SA)))==FALSE)
RMcont<- which(is.na(match(colnames(rfimps),colnames(RMvars)))==FALSE)
Clincont<-which(is.na(match(colnames(rfimps),colnames(Clin)))==FALSE)

S<- (dim(peaks)[1]*3)
classimp=matrix(0,S,8)
for(i in 1:S){
	classimp[i,1]<- mean(as.numeric(rfimps[i,CNVcont]),na.rm=TRUE)
	classimp[i,2]<- mean(as.numeric(rfimps[i,Somacont]),na.rm=TRUE)
	classimp[i,3]<- mean(as.numeric(rfimps[i,TFcont]),na.rm=TRUE)
	classimp[i,4]<- mean(as.numeric(rfimps[i,Mcont]),na.rm=TRUE)
	classimp[i,5]<- mean(as.numeric(rfimps[i,SGcont]),na.rm=TRUE)
	classimp[i,6]<- mean(as.numeric(rfimps[i,SAcont]),na.rm=TRUE)
	classimp[i,7]<- mean(as.numeric(rfimps[i,RMcont]),na.rm=TRUE)
	classimp[i,8]<- mean(as.numeric(rfimps[i,Clincont]),na.rm=TRUE)
	
	}

colnames(classimp)<-c("CNV","MUT","TF","MEC","SGOC","SAM","RM","CLINICAL")
head(classimp)

write.table(classimp
,"average var imp from classes lgg  low mse modeling in lgg.txt")

elcoefs<- nco[gc,-1]

CNVcon<- which(is.na(match(colnames(elcoefs),colnames(cnvMat)))==FALSE)
Somacon<-which(is.na(match(colnames(elcoefs),colnames(MutMat)))==FALSE)
TFcon<-which(is.na(match(colnames(elcoefs),colnames(TFvars)))==FALSE)
Mcon<-which(is.na(match(colnames(elcoefs),colnames(MEC)))==FALSE)
SGcon<-which(is.na(match(colnames(elcoefs),colnames(SGOCvars)))==FALSE)
SAcon<-which(is.na(match(colnames(elcoefs),colnames(SA)))==FALSE)
RMcon<- which(is.na(match(colnames(elcoefs),colnames(RMvars)))==FALSE)
Clincon<-which(is.na(match(colnames(elcoefs),colnames(Clin)))==FALSE)


classcoef=matrix(0,S,8)
for(i in 1:S){
	classcoef[i,1]<- mean(as.numeric(abs(elcoefs[i,CNVcon])),na.rm=TRUE)
	classcoef[i,2]<- mean(as.numeric(abs(elcoefs[i,Somacon])),na.rm=TRUE)
	classcoef[i,3]<- mean(as.numeric(abs(elcoefs[i,TFcon])),na.rm=TRUE)
	classcoef[i,4]<- mean(as.numeric(abs(elcoefs[i,Mcon])),na.rm=TRUE)
	classcoef[i,5]<- mean(as.numeric(abs(elcoefs[i,SGcon])),na.rm=TRUE)
	classcoef[i,6]<- mean(as.numeric(abs(elcoefs[i,SAcon])),na.rm=TRUE)
	classcoef[i,7]<- mean(as.numeric(abs(elcoefs[i,RMcon])),na.rm=TRUE)
	classcoef[i,8]<- mean(as.numeric(abs(elcoefs[i,Clincon])),na.rm=TRUE)
	
	}

colnames(classcoef)<-c("CNV","MUT","TF","MEC","SGOC","SAM","RM","CLINICAL")
head(classcoef)

write.table(classcoef
,"fraction variable usage from classes lgg low mse modeling in lgg.txt")

Elsum=matrix(0,8,3)

Elsum[,1]<- (colSums(classcoef)/S)*100

for(i in 1:8){
	Elsum[i,2]<- max(classcoef[,i])*100
	Elsum[i,3]<- min(classcoef[,i])*100
}


RFsum=matrix(0,8,3)

RFsum[,1]<- colSums(classimp)/S

for(i in 1:8){
	RFsum[i,2]<- max(classimp[,i])
	RFsum[i,3]<- min(classimp[,i])
}

write.table(Elsum,"summary elastic net classes in lgg.txt")
write.table(RFsum,"summary random forest classes in lgg.txt")

# ********* #
#  FIGURE 3 #
# ********* #


 # Pick best MEC gene based on correlation in fig 1 and 2


# Spearman Cors for all probes with the top MEC gene (BHMT2 in LGG):

BHMT2<- nLGGEXP[which(nGenenames=="BHMT2"),]

dim(TAB)

cor=c()
pval=c()
for (i in 1:(dim(TAB)[1])){
	cor[i]<- cor.test(BHMT2,as.numeric(TAB[i,12:(N+11)]), method="spearman")$estimate
	pval[i]<- cor.test(BHMT2,as.numeric(TAB[i,12:(N+11)]), method="spearman")$p.value
	cat(i)
}

setwd("/home/X/LGG")
write.table(cbind(cor,pval),"all sorted probes BHMT2 cors in LGG.txt")

# sliding window

setwd("/home/X/LGG")
corspvals<- read.table("all sorted probes BHMT2 cors in LGG.txt")
dim(corspvals)
cor<- as.numeric(corspvals[,1])
pval<- as.numeric(corspvals[,2])

chroms<- as.numeric(unique(TAB$chr))
speaks<- matrix(0,1,3)
colnames(speaks)<- c("chr","start","end")


for(i in 1:length(chroms)) {
	ind<- which(TAB$chr==chroms[i])
	tab<- cbind(TAB[ind,1:11], cor[ind], pval[ind])
	colnames(tab)[12:13]<- c("cor","pval")
	thresh<- as.numeric(quantile(abs(cor[ind]),0.90))
# find the top 10% probes with highest correlation across the genome	
	hits<- which(abs(tab$cor)>thresh)
	speaks<- rbind(speaks,c(chroms[i],0,0))
	 for(j in 1:length(hits)){
	 	location<- tab$position[hits[j]]
	# analyze all probes in a 6kb window surrounding the top probes
	 	if(location>speaks[dim(speaks)[1],3]){
	 		window<- which(abs(tab$position-location)<3000)
	 		good<- which(tab$cor[window]*tab$cor[hits[j]]>0)
	 		if(length(good)>2){
	 		finalgood<- length(which(tab$pval[window[good]]<0.00001))
	 		frac<- finalgood/length(window)
	 
	 # select windows in which 80% of the probes have a significant correlation
	 		if(frac >= 0.8){
	 			newpeak<- c(as.numeric(chroms[i]), as.numeric(tab$position[window[1]]), as.numeric(tab$position[window[length(window)]]))
	 			speaks<- rbind(speaks,newpeak)
	 			}
	 		}
	 	}
	 	
	 	
	 }
	 cat(i)
}

write.table(speaks, "all BHMT2 peaks in LGG using sliding window.txt")

dim(speaks)
    

# annotate identified regions

GAP=matrix(0,1,13)
colnames(GAP)[12:13]<- c("cor","pval")
colnames(GAP)[1:11]<- colnames(TAB)[1:11]


peakannot=matrix(0,1,13)
colnames(peakannot)[12:13]<- c("cor","pval")
colnames(peakannot)[1:11]<- colnames(TAB)[1:11]

for(i in 1:(dim(speaks)[1])){
	if(speaks[i,3]>0){
		CI<- which(TAB$chr==as.character(speaks[i,1]))
		nTAB<- TAB[CI,1:11]
		nBHMT2<- corspvals[CI,]
		start<- which(nTAB$position==speaks[i,2])
		end<- which(nTAB$position==speaks[i,3])
		pk<- cbind(nTAB[start:end,],nBHMT2[start:end,])
		peakannot<- rbind(peakannot,pk)
		peakannot<- rbind(peakannot,GAP)
			}
	
}

table(peakannot$chr)

peakannot[1:100,1:5]


setwd("/home/X/LGG")
write.table(peakannot, "peak annotatioins all peaks from sliging window for bhmt2 in lgg.txt")


# pathway enrichment analyses was done using EnrichR

# distance to nearest TSS was found using GREAT


 # fig S4 : test of Met cycle specificity:


speaks<- read.table("all BHMT2 peaks in LGG using sliding window.txt", row.names=NULL)
speaks<- speaks[,-1]
peakannot<- read.table("peak annotatioins all peaks from sliging window for bhmt2 in lgg.txt")

# select top regions that form "peaks" based on information theory peak shape

library("pastecs")

score=c()
medcor=c()
for(i in 1:(dim(speaks)[1])){
	if(speaks[i,2]>0){
		ch<- which(peakannot$chr==speaks[i,1])
		tab<- peakannot[ch,]
		start<- which(tab$position==speaks[i,2])
		end<- which(tab$position==speaks[i,3])
		corvec<- tab$cor[start:end]
		medcor[i]<- median(abs(corvec))
		PeakPit<- turnpoints(corvec)
        score[i] <- PeakPit$proba
	}
	
}


top<- which(medcor>0.4)
shape<- which(score[top]<0.1)

hits<- speaks[top[shape],]
 write.table(hits,"pastec 141 good peaks for bhmt2 in lgg.txt")


  #calculate window test randomization p-val for each of the top peaks 

BHMT2<- nLGGEXP[which(nGenenames=="BHMT2"),]

hits<- as.data.frame(hits)

randpval=c()
for(i in 1:(dim(hits)[1])){
ch<- hits$chr[i]
tab<- TAB[which(TAB$chr==ch),]
start<- which(tab$position==hits$start[i])
end<- which(tab$position==hits$end[i])
final<- tab[start:end,]
f<- dim(final)[1]
BHcors=c()
for(n in 1:f){
	BHcors[n]<- cor.test(BHMT2,as.numeric(final[n,12:(N+11)]), method="spearman")$estimate
}
BH<- length(which(abs(BHcors)>0.3))
# correlate with 500 random genes
random<- sample(seq(1,dim(nLGGEXP)[1]),500)
randcors=c()
 for(j in 1:500){
 	rcors=c()
for(n in 1:f){
	rcors[n]<- cor.test(nLGGEXP[random[j],],as.numeric(final[n,12:(N+11)]), method="spearman")$estimate
}
randcors[j]<- length(which(abs(rcors)>0.3))

 }
randpval[i]<- length(which(randcors>BH))/500

cat(i)
}

randpval

setwd("/home/X/LGG")
write.table(randpval,"lgg randomization p-values for 141 bhmt2 peaks in the window analysis.txt")

# ******** #
# FIGURE 4 # 
# ******** #


cancerdrivers<- read.table("NOV cancer gene list.txt")
cancergenes<- as.character(cancerdrivers[,1])


tind<- which(match(table$UCSCRefGeneName, cancergenes)>0)

GTAB<- cbind(allanot[tind,],table[tind,],METHYLATION[tind,])


GTAB<- GTAB[with(GTAB,order(chr,position)),]

Q<- unique(GTAB$UCSCRefGeneName)
N<- dim(METHYLATION)[2]

# calculated average gene body and promoter methylation for each of the cancer genes 
cancerResponse=matrix(0,1,N+2)

for(i in 1:length(Q)){
	name<- which(GTAB$UCSCRefGeneName==Q[i])
	new<- GTAB[name,]
	promoter<- c(which(new$UCSCRefGeneGroup=="TSS1500"),which(new$UCSCRefGeneGroup=="TSS200"), which(new$UCSCRefGeneGroup=="5'UTR"))
	body<- c(which(new$UCSCRefGeneGroup=="Body"),which(new$UCSCRefGeneGroup=="3'UTR"),which(new$UCSCRefGeneGroup=="1stExon"))
	
	meanP=c()
	meanB=c()
	for(j in 1:N){
		meanP[j]<- mean(new[promoter,j+11],na.rm=TRUE)
		meanB[j]<- mean(new[body,j+11],na.rm=TRUE)
	}
	
P<- c(as.character(Q[i]),"promoter",meanP)
B<- c(as.character(Q[i]),"body",meanB)
cancerResponse<- rbind(cancerResponse,P)
cancerResponse<- rbind(cancerResponse,B)	
cat(i)	
}


cancerResponse<- cancerResponse[-1,]

setwd("/home/X/LGG")
write.table(cancerResponse,"FINAL cancer genes methyaltion in lgg.txt")

  # modeling:# just using MeC variables :
  
library(randomForest)
library(glmnet)  

cancerResponse<- read.table("FINAL cancer genes methyaltion in lgg.txt", row.names=NULL) 
cancerResponse<- cancerResponse[,-1]
  
df.feat<- data.frame(t(nLGGEXP[MECEXP,]))
colnames(df.feat)<- as.character(nGenenames[MECEXP])

# exclude NA rows and extra columns

response<- data.matrix(cancerResponse[-34,-c(1,2)])

# set cross validation fold-k

k=7

test.ind.mat <- matrix(sample(1:dim(response)[2]),nrow=k)

rfMSE.mat <- matrix(0,nrow=dim(response)[1],ncol=k)
elMSE.mat <- matrix(0,nrow=dim(response)[1],ncol=k)
varimp<- matrix(0,1,dim(df.feat)[2])
varcoef<- matrix(0,1,(dim(df.feat)[2]+1))
rfpredicted<- matrix(0,dim(response)[1],N)
elpredicted<- matrix(0,dim(response)[1],N)

for (i in 1:k) {
	
	# Run random forest and elastic net on each response
	for (j in 1:dim(response)[1]) {
		
		df.resp <- response[j,]
df.full <- data.frame("resp"=as.numeric(df.resp),df.feat)

  # Create train and test set
	df.train <- df.full[-test.ind.mat[i,],]
	df.test <- df.full[test.ind.mat[i,],]
	
	# use na.rough to impute missing values
	df.train.rough <- na.roughfix(df.train)
	df.test.rough <- na.roughfix(df.test)
	

   #rf
		rf <- randomForest(y=df.train.rough[,1],x=df.train.rough[,-1],ntree=500,mtry=3,importance=TRUE)
		rf.pred <- predict(rf,df.test.rough)
		rfMSE.mat[j,i] <- (1/(dim(test.ind.mat)[2]))*(sum((rf.pred - df.test.rough[,1])^2))
		varimp<- rbind(varimp,importance(rf,scale=TRUE)[,1])
rfpredicted[j,test.ind.mat[i,]]<- rf.pred


    #el
    

Y<- df.train.rough[,1]

Ecvfit = cv.glmnet(data.matrix(df.train.rough[,-1]),Y,alpha=0.5,nfold=5)


tfit<- predict(Ecvfit,newx=data.matrix(df.test.rough[,-1]) ,s="lambda.min")

elMSE.mat[j,i] <- (1/(dim(test.ind.mat)[2]))*(sum((tfit - df.test.rough[,1])^2))
varcoef<- rbind(varcoef,coef(Ecvfit, s = "lambda.min")[,1])
elpredicted[j,test.ind.mat[i,]]<-tfit

				
		cat("Fold:  ",i,"  of 7;               response ",j," . \n")
	}
	
}

setwd("/home/X/LGG")
write.table(varimp,"final cancer genes MEC only var imps in lgg.txt")
write.table(varcoef,"final cancer genes MEC only var coefs in lgg.txt")
write.table(rfMSE.mat,"final cancer genes MEC only RF MSE mat in lgg.txt")
write.table(elMSE.mat,"final cancer genes MEC only EL MSE mat in lgg.txt")
write.table(rfpredicted,"final cancer genes MEC only RF prediction matrix in lgg.txt")
write.table(elpredicted,"final cancer genes MEC only EL prediction matrix in lgg.txt")


mean(data.matrix(rfMSE.mat))

mean(data.matrix(elMSE.mat))


rfE<- sqrt(rfMSE.mat)
elE<- sqrt(elMSE.mat)

rfError<- rowSums(rfE)/k
elError<- rowSums(elE)/k


length(which(rfError<0.1))

length(which(elError<0.1))
  


vim<- varimp[-1,]
meanImp<- colSums (vim)/(dim(vim)[1])

meanImp
   
vco<- varcoef[-1,]
meanCo<- colSums (vco)/(dim(vim)[1])

meanCo
nco<- data.matrix(vco)
nco[which(nco!=0)]<- 1
varuse<- colSums(nco)

varuse/(k*dim(response)[1])


# COMPARE CANCER GENES WITH PERMUTATED METHYLATIONS OR RANDOM NUMBERS

 #permutations:
 
scramble <- function(vec,n) {
  res <- tapply(vec,(seq_along(vec)+n-1)%/%n,
                FUN=function(x) x[sample.int(length(x), size=length(x))])
  unname(unlist(res))
}

N
z=dim(response)[1]
randomResponse=matrix(0,z,N)

for(i in 1:z){
	randomResponse[i,]<- scramble(response[i,],N)
	
}



k

RrfMSE.mat <- matrix(0,nrow=dim(randomResponse)[1],ncol=k)
RelMSE.mat <- matrix(0,nrow=dim(randomResponse)[1],ncol=k)
Rvarimp<- matrix(0,1,dim(df.feat)[2])
Rvarcoef<- matrix(0,1,(dim(df.feat)[2]+1))

for (i in 1:k) {
	
	# Run random forest and elastic net on each response
	for (j in 1:dim(randomResponse)[1]) {
		
		df.resp <- randomResponse[j,]
df.full <- data.frame("resp"=as.numeric(df.resp),df.feat)

  # Create train and test set
	df.train <- df.full[-test.ind.mat[i,],]
	df.test <- df.full[test.ind.mat[i,],]
	
	# use na.rough to impute missing values
	df.train.rough <- na.roughfix(df.train)
	df.test.rough <- na.roughfix(df.test)
	

   #rf
		rf <- randomForest(y=df.train.rough[,1],x=df.train.rough[,-1],ntree=500,mtry=3,importance=TRUE)
		rf.pred <- predict(rf,df.test.rough)
		RrfMSE.mat[j,i] <- (1/(dim(test.ind.mat)[2]))*(sum((rf.pred - df.test.rough[,1])^2))
		Rvarimp<- rbind(Rvarimp,importance(rf,scale=TRUE)[,1])

    #el
    

Y<- df.train.rough[,1]

Ecvfit = cv.glmnet(data.matrix(df.train.rough[,-1]),Y,alpha=0.5,nfold=5)


tfit<- predict(Ecvfit,newx=data.matrix(df.test.rough[,-1]) ,s="lambda.min")

RelMSE.mat[j,i] <- (1/(dim(test.ind.mat)[2]))*(sum((tfit - df.test.rough[,1])^2))
Rvarcoef<- rbind(Rvarcoef,coef(Ecvfit, s = "lambda.min")[,1])
				
		cat("Fold:  ",i,"  of 7;               response ",j," . \n")
	}
	
}

setwd("/home/X/LGG") 
write.table(RrfMSE.mat,"permutation rf mse from MEC only in lgg.txt")
write.table(RelMSE.mat,"permutation el mse from MEC only in lgg.txt")
write.table(Rvarimp,"permutation rf var imps from MEC only in lgg.txt")
write.table(Rvarcoef,"permutation el var coefs from MEC only in lgg.txt")



vim<- Rvarimp[-1,]
meanImp<- colSums (vim)/(dim(vim)[1])
meanImp
     

vco<- Rvarcoef[-1,]

nco<- vco
nco[which(nco!=0)]<- 1
varuse<- colSums(nco)

varuse/(45*k)



ks.test(as.vector(varimp),as.vector(Rvarimp))



ks.test(as.vector(varcoef[,-1]),as.vector(Rvarcoef[,-1]))


   # Random numbers:
   
  z=dim(response)[1]
randomResponse=matrix(0,z,N)

for(i in 1:z){
	randomResponse[i,]<- runif(N,0,1)
	
}

k

RrfMSE.mat <- matrix(0,nrow=dim(randomResponse)[1],ncol=k)
RelMSE.mat <- matrix(0,nrow=dim(randomResponse)[1],ncol=k)
Rvarimp<- matrix(0,1,dim(df.feat)[2])
Rvarcoef<- matrix(0,1,(dim(df.feat)[2]+1))

for (i in 1:k) {
	
	# Run random forest and elastic net on each response
	for (j in 1:dim(randomResponse)[1]) {
		
		df.resp <- randomResponse[j,]
df.full <- data.frame("resp"=as.numeric(df.resp),df.feat)

  # Create train and test set
	df.train <- df.full[-test.ind.mat[i,],]
	df.test <- df.full[test.ind.mat[i,],]
	
	# use na.rough to impute missing values
	df.train.rough <- na.roughfix(df.train)
	df.test.rough <- na.roughfix(df.test)
	

   #rf
		rf <- randomForest(y=df.train.rough[,1],x=df.train.rough[,-1],ntree=500,mtry=3,importance=TRUE)
		rf.pred <- predict(rf,df.test.rough)
		RrfMSE.mat[j,i] <- (1/(dim(test.ind.mat)[2]))*(sum((rf.pred - df.test.rough[,1])^2))
		Rvarimp<- rbind(Rvarimp,importance(rf,scale=TRUE)[,1])

    #el
    

Y<- df.train.rough[,1]

Ecvfit = cv.glmnet(data.matrix(df.train.rough[,-1]),Y,alpha=0.5,nfold=5)


tfit<- predict(Ecvfit,newx=data.matrix(df.test.rough[,-1]) ,s="lambda.min")

RelMSE.mat[j,i] <- (1/(dim(test.ind.mat)[2]))*(sum((tfit - df.test.rough[,1])^2))
Rvarcoef<- rbind(Rvarcoef,coef(Ecvfit, s = "lambda.min")[,1])
				
		cat("Fold:  ",i,"  of 7;               response ",j," . \n")
	}
	
}


setwd("/home/X/LGG") 
write.table(RrfMSE.mat,"random numbers rf mse from MEC only in lgg.txt")
write.table(RelMSE.mat,"random numbers el mse from MEC only in lgg.txt")
write.table(Rvarimp,"random numbers rf var imps from MEC only in lgg.txt")
write.table(Rvarcoef,"random numbers el var coefs from MEC only in lgg.txt")

mean(RrfMSE.mat)
mean(RelMSE.mat)
 

rfE<- sqrt(RrfMSE.mat)
elE<- sqrt(RelMSE.mat)

rfError<- rowSums(rfE)/k
elError<- rowSums(elE)/k


vim<- Rvarimp[-1,]
meanImp<- colSums (vim)/(dim(vim)[1])

meanImp
     
      
vco<- Rvarcoef[-1,]

nco<- vco
nco[which(nco!=0)]<- 1
varuse<- colSums(nco)

varuse/(45*k)

ks.test(as.vector(rfMSE.mat),as.vector(RrfMSE.mat))

ks.test(as.vector(elMSE.mat),as.vector(RelMSE.mat))
 
  # mse comparisons
  
randnoise=read.table("random numbers el mse from MEC only in lgg.txt")
permut=read.table("permutation el mse from MEC only in lgg.txt")
cancer=read.table("final cancer genes MEC only EL MSE mat in lgg.txt")


write.table(cbind(as.vector(unlist(cancer)), as.vector(unlist(permut)),as.vector(unlist(randnoise))), "lgg mse of mec comparison between random sets in prism.txt")


# variable usage comparisons
randnoise=read.table("random numbers el var coefs from MEC only in lgg.txt")
permut=read.table("permutation el var coefs from MEC only in lgg.txt")
cancer=read.table("final cancer genes MEC only var coefs in lgg.txt")


randnoise<- data.matrix(randnoise[-1,-1])
randnoise[which(randnoise!=0)]<- 1
permut<- data.matrix(permut[-1,-1])
permut[which(permut!=0)]<- 1
cancer<- data.matrix(cancer[-1,-1])
cancer[which(cancer!=0)]<- 1
write.table(cbind(as.vector(unlist(cancer)), as.vector(unlist(permut)),as.vector(unlist(randnoise))), "lgg var usage of mec comparison between random sets in prism.txt")


 
# MODELING CANCER GENES USING ALL VARIABLES:

setwd("/home/X/LGG")
cancerResponse<- read.table("FINAL cancer genes methyaltion in lgg.txt", row.names=NULL) 
cancerResponse<- cancerResponse[,-1]
response<- cancerResponse[,-c(1,2)]

Clin<- read.csv("small clinical variables for lgg.csv")
Clin<- Clin[,-1]
SAMvars<- read.table("SAM metabolizing variables in lgg.txt")
RMvars<- read.table("new chromatin remodeller and dna meths variables in lgg.txt")
SGOCvars<- read.table("SGOC variables in lgg.txt")
MEC<-SGOCvars[,c(4,6,8,9)]
SGOCvars<- SGOCvars[,-c(4,6,8,9)]

TFvars<- read.table("small list transcription factor variables in lgg.txt")
MutMat<- read.table("final mutations matrix in lgg.txt")
cnvMat<- read.table("final cnvs matrix in lgg.txt")


X<- cbind(Clin,MutMat,cnvMat,MEC, TFvars,RMvars,SAMvars,SGOCvars)


X1<- X[, unique(colnames(X))]
dim(X1)

W<- dim(X1)[2]

# set cross validation folds

k=3

test.ind.mat <- matrix(sample(1:dim(response)[2]),nrow=k)

rfMSE.mat <- matrix(0,nrow=dim(response)[1],ncol=k)
elMSE.mat <- matrix(0,nrow=dim(response)[1],ncol=k)
varimp<- matrix(0,1,W)
varcoef<- matrix(0,1,W+1)
rfpredicted<- matrix(0,dim(response)[1],N)
elpredicted<- matrix(0,dim(response)[1],N)



for (i in 1:k) {
	
	# Run random forest and elastic net on each response
	for (j in 1:dim(response)[1]) {
		
		df.resp <- response[j,]
df.full <- data.frame("resp"=as.numeric(df.resp),X1)

  # Create train and test set
	df.train <- df.full[-test.ind.mat[i,],]
	df.test <- df.full[test.ind.mat[i,],]
	
	# use na.rough to impute missing values
	df.train.rough <- na.roughfix(df.train)
	df.test.rough <- na.roughfix(df.test)
	

   #rf
		rf <- randomForest(y=df.train.rough[,1],x=df.train.rough[,-1],ntree=500,importance=TRUE)
		rf.pred <- predict(rf,df.test.rough)
		rfMSE.mat[j,i] <- (1/(dim(test.ind.mat)[2]))*(sum((rf.pred - df.test.rough[,1])^2))
		varimp<- rbind(varimp,importance(rf,scale=TRUE)[,1])
rfpredicted[j,test.ind.mat[i,]]<- rf.pred

    #el
    

Y<- df.train.rough[,1]

Ecvfit = cv.glmnet(data.matrix(df.train.rough[,-1]),Y,alpha=0.5,nfold=5)


tfit<- predict(Ecvfit,newx=data.matrix(df.test.rough[,-1]) ,s="lambda.min")

elMSE.mat[j,i] <- (1/(dim(test.ind.mat)[2]))*(sum((tfit - df.test.rough[,1])^2))
varcoef<- rbind(varcoef,coef(Ecvfit, s = "lambda.min")[,1])
		elpredicted[j,test.ind.mat[i,]]<-tfit		
		cat("Fold:  ",i,"  of 3;               response ",j," . \n")
	}
	
}

setwd("/home/X/LGG")
write.table(varimp,"final cancer genes ALL VARS var imps in lgg.txt")
write.table(varcoef,"final cancer genes ALL VARS var coefs in lgg.txt")
write.table(rfMSE.mat,"final cancer genes ALL VARS RF MSE mat in lgg.txt")
write.table(elMSE.mat,"final cancer genes ALL VARS EL MSE mat in lgg.txt")
write.table(rfpredicted,"final cancer genes ALL VARS RF prediction matrix in lgg.txt")
write.table(elpredicted,"final cancer genes ALL VARS EL prediction matrix in lgg.txt")

# prediction errors

mean(data.matrix(elMSE.mat),na.rm=TRUE)
mean(data.matrix(rfMSE.mat),na.rm=TRUE)

rfE<- sqrt(rfMSE.mat)
elE<- sqrt(elMSE.mat)

rfError<- rowSums(rfE)/k
elError<- rowSums(elE)/k

length(which(rfError<0.1))

length(which(elError<0.1))
 


vim<- varimp
meanImp<- colSums (vim)/(length(names)*k)

sort(meanImp)
     
write.table(meanImp,"cancer genes mean variable importances all variables models in lgg.txt") 


    
vco<- varcoef
nco<- data.matrix(vco)

nco[which(nco!=0)]<- 1
varuse<- colSums(nco)

sort(varuse)

write.table(varuse,"cancer genes variable usage all variables models in lgg.txt")


# MEC vars contributions:

MeC<- c(29,30,31,32)
meanImp[MeC]

rank(meanImp)[MeC]/W
 

varuse[MeC]/((length(names))*k)

     
Gvar<- varuse

rank(Gvar)[MeC]/W

  
# Contribution of variable CLASSES (average)

which(colnames(SAMvars)=="MAT2B") 
SA<- SAMvars[,-40] 

rfimps<- varimp

CNVcont<- which(is.na(match(colnames(rfimps),colnames(cnvMat)))==FALSE)
Somacont<-which(is.na(match(colnames(rfimps),colnames(MutMat)))==FALSE)
TFcont<-which(is.na(match(colnames(rfimps),colnames(TFvars)))==FALSE)
Mcont<-which(is.na(match(colnames(rfimps),colnames(MEC)))==FALSE)
SGcont<-which(is.na(match(colnames(rfimps),colnames(SGOCvars)))==FALSE)
SAcont<-which(is.na(match(colnames(rfimps),colnames(SA)))==FALSE)
RMcont<- which(is.na(match(colnames(rfimps),colnames(RMvars)))==FALSE)
Clincont<-which(is.na(match(colnames(rfimps),colnames(Clin)))==FALSE)


classimp=matrix(0,(dim(rfimps)[1]),8)
for(i in 1:(dim(rfimps)[1])){
	classimp[i,1]<- mean(as.numeric(rfimps[i,CNVcont]),na.rm=TRUE)
	classimp[i,2]<- mean(as.numeric(rfimps[i,Somacont]),na.rm=TRUE)
	classimp[i,3]<- mean(as.numeric(rfimps[i,TFcont]),na.rm=TRUE)
	classimp[i,4]<- mean(as.numeric(rfimps[i,Mcont]),na.rm=TRUE)
	classimp[i,5]<- mean(as.numeric(rfimps[i,SGcont]),na.rm=TRUE)
	classimp[i,6]<- mean(as.numeric(rfimps[i,SAcont]),na.rm=TRUE)
	classimp[i,7]<- mean(as.numeric(rfimps[i,RMcont]),na.rm=TRUE)
	classimp[i,8]<- mean(as.numeric(rfimps[i,Clincont]),na.rm=TRUE)
	
	}


write.table(classimp
,"cancer genes average var imp from classes lgg 3*45 modeling.txt")

elcoefs<- nco

CNVcon<- which(is.na(match(colnames(elcoefs),colnames(cnvMat)))==FALSE)
Somacon<-which(is.na(match(colnames(elcoefs),colnames(MutMat)))==FALSE)
TFcon<-which(is.na(match(colnames(elcoefs),colnames(TFvars)))==FALSE)
Mcon<-which(is.na(match(colnames(elcoefs),colnames(MEC)))==FALSE)
SGcon<-which(is.na(match(colnames(elcoefs),colnames(SGOCvars)))==FALSE)
SAcon<-which(is.na(match(colnames(elcoefs),colnames(SA)))==FALSE)
RMcon<- which(is.na(match(colnames(elcoefs),colnames(RMvars)))==FALSE)
Clincon<-which(is.na(match(colnames(elcoefs),colnames(Clin)))==FALSE)


classcoef=matrix(0,(dim(elcoefs)[1]),8)
for(i in 1:(dim(elcoefs)[1])){
	classcoef[i,1]<- mean(as.numeric(abs(elcoefs[i,CNVcon])),na.rm=TRUE)
	classcoef[i,2]<- mean(as.numeric(abs(elcoefs[i,Somacon])),na.rm=TRUE)
	classcoef[i,3]<- mean(as.numeric(abs(elcoefs[i,TFcon])),na.rm=TRUE)
	classcoef[i,4]<- mean(as.numeric(abs(elcoefs[i,Mcon])),na.rm=TRUE)
	classcoef[i,5]<- mean(as.numeric(abs(elcoefs[i,SGcon])),na.rm=TRUE)
	classcoef[i,6]<- mean(as.numeric(abs(elcoefs[i,SAcon])),na.rm=TRUE)
	classcoef[i,7]<- mean(as.numeric(abs(elcoefs[i,RMcon])),na.rm=TRUE)
	classcoef[i,8]<- mean(as.numeric(abs(elcoefs[i,Clincon])),na.rm=TRUE)
	
	}



write.table(classcoef
,"cancer genes fraction variable usage from classes lgg 3*45 modeling.txt")


# collapse class contributions

classRF<- matrix(0,45,8)
classEL<- matrix(0,45,8)
for(i in 1:45){
	
	classRF[i,]<- (classimp[i,]+classimp[i+45,]+classimp[i+90,])/3
	classEL[i,]<- (classcoef[i,]+classcoef[i+45,]+classcoef[i+90,])/3
}
colnames(classRF)<- c("CNV","MUT","TF","MEC","SGOC","SAM","RM","CLINICAL")

colnames(classEL)<- c("CNV","MUT","TF","MEC","SGOC","SAM","RM","CLINICAL")

rownames(classRF)<- names
rownames(classEL)<- names

write.table(classRF,"average class contributions in lgg from RF in cancer genes for heatmaps.txt")

write.table(classEL,"average class contributions in lgg from EL in cancer genes for heatmaps.txt")


 Elsum=matrix(0,8,3)

Elsum[,1]<- (colSums(classcoef)/(3*45))*100

for(i in 1:8){
	Elsum[i,2]<- max(classcoef[,i])*100
	Elsum[i,3]<- min(classcoef[,i])*100
}


RFsum=matrix(0,8,3)

RFsum[,1]<- colSums(classimp)/(3*45)

for(i in 1:8){
	RFsum[i,2]<- max(classimp[,i])
	RFsum[i,3]<- min(classimp[,i])
}

write.table(Elsum,"cancer genes summary elastic net classes in lgg.txt")
write.table(RFsum,"cancer genes summary random forest classes in lgg.txt")

ELlgg<- read.table("final cancer genes ALL VARS var coefs in lgg.txt")
ELlgg<- data.matrix(ELlgg[-1,])
ELlgg[which(ELlgg!=0)]<- 1
lggvaruse<- colSums(ELlgg)
write.table(sort(lggvaruse), "STAR lgg ranks for prism.txt")

# ******** #
# Figure 5 #    
# ******** #


# download survival data using tcga assembler:

install.packages(c("HGNChelper", "RCurl", "httr", "stringr", "digest", "bitops"), dependencies=T) 

setwd("/home/X/TCGAAssembler") 

source("./Module_A.r"); 
source("./Module_B.r"); 

lggclin<- DownloadClinicalData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = "./lggsurvival", cancerType = "LGG", clinicalDataType = c("patient", "drug", "days_to_death" ,"days_to_last_follow_up", "follow_up")); 


setwd("/home/X/TCGAAssembler/lggsurvival")
list.files("/home/X/TCGAAssembler/lggsurvival")
follow<- read.table("nationwidechildrens.org_clinical_follow_up_v1.0_lgg.txt",sep="\t")
 dim(follow)
 
 follow<- follow[-1,]
 colnames(follow)<- as.character(unlist(follow[1,]))
 
 follow<- follow[-1,]
 colnames(follow)
 
 follow$days_to_death
 follow$vital_status
 

  # make recurrance matrix:
  
bars<- follow$bcr_patient_barcode

length(bars)
length(unique(bars))

setwd("/home/X/LGG")
 barcodes<- read.table("LGG BARCODES 534.txt")

sinds=c()
for(i in 1:N){
 	if(is.na(match(substr(barcodes[i,1],1,12),bars))==FALSE){
 		
 	ind<- which(bars==substr(barcodes[i,1],1,12))
 	
 	if(length(ind)>1){
 		
 		 			ind<- ind[1]
 		}
 		sinds[i]<- ind
 		}}

SurMat=follow[sinds[1],]
blank=matrix(0,1,dim(follow)[2])
blank<- as.data.frame(blank)
colnames(blank)<- colnames(follow)

for(i in 2:N){
	if(is.na(sinds[i])==FALSE){
		SurMat<- rbind(SurMat, follow[sinds[i],])
	}
	else {
		SurMat<- rbind(SurMat,blank)
	}
}



setwd("/home/X/LGG")
write.table(SurMat,"survival matrix in lgg.txt")
   
  
# use death

SurMat<- read.table("survival matrix in lgg.txt")
SUR<- SurMat$days_to_death
alive<- grep("[Not Applicable]",SurMat$days_to_death)
SUR<- as.vector(SUR)


cancerResponse<- read.table("FINAL cancer genes methyaltion in lgg.txt", row.names=NULL) 
cancerResponse<- cancerResponse[,-1]
response<- cancerResponse[,-c(1,2)]


library("survival")

# separate samples based on the prediction errors of models of cancer genes only using Met cycle


elpredicted<- read.table("final cancer genes MEC only EL prediction matrix in lgg.txt")

dif<- abs(response-elpredicted)
sampPred<- colSums(dif)/45
range(sampPred)
hist(sampPred)
med=0.04

sampLab=c()
sampLab[which(sampPred<med)]<- "good"
sampLab[which(sampPred>med)]<- "bad"

test<- data.frame("time"=as.numeric(SUR),"status"=sampLab)

# CENSOR DATA:

clinical<- SurMat

ind_keep <- grep("days_to_new_tumor_event_after_initial_treatment",colnames(clinical))

# recurrence
new_tum <- as.matrix(clinical[,ind_keep])
new_tum_collapsed <- c()
for (i in 1:dim(new_tum)[1]){
  if(sum(is.na(new_tum[i,])) < dim(new_tum)[2]){
    m <- max(new_tum[i,],na.rm=T)
    new_tum_collapsed <- c(new_tum_collapsed,m)
  } else {
    new_tum_collapsed <- c(new_tum_collapsed,"NA")
  }
}

# do the same to death
ind_keep <- grep("days_to_death",colnames(clinical))
death <- as.matrix(clinical[,ind_keep])
death_collapsed <- c()
for (i in 1:dim(death)[1]){
  if(sum(is.na(death[i,])) < dim(death)[2]){
    m <- max(death[i,],na.rm=T)
    death_collapsed <- c(death_collapsed,m)
  } else {
    death_collapsed <- c(death_collapsed,"NA")
  }
}

# and days last follow up here we take the most recent which is the min number
ind_keep <- grep("days_to_last_followup",colnames(clinical))
fl <- as.matrix(clinical[,ind_keep])
fl_collapsed <- c()
for (i in 1:dim(fl)[1]){
  if(sum(is.na(fl[i,])) < dim(fl)[2]){
    m <- min(fl[i,],na.rm=T)
    fl_collapsed <- c(fl_collapsed,m)
  } else {
    fl_collapsed <- c(fl_collapsed,"NA")
  }
}

# and put everything together
all_clin <- data.frame(new_tum_collapsed,death_collapsed,fl_collapsed)
colnames(all_clin) <- c("new_tumor_days", "death_days", "followUp_days")


all_clin$new_time <- c()
for (i in 1:length(as.numeric(as.character(all_clin$new_tumor_days)))){
  all_clin$new_time[i] <- ifelse(is.na(as.numeric(as.character(all_clin$new_tumor_days))[i]),
                    as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$new_tumor_days))[i])
}

# create vector time to death containing values to censor for death

all_clin$new_death <- c()
for (i in 1:length(as.numeric(as.character(all_clin$death_days)))){
  all_clin$new_death[i] <- ifelse(is.na(as.numeric(as.character(all_clin$death_days))[i]),
                                 as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$death_days))[i])
}

# create vector for death censoring
table(clinical$vital_status)
all_clin$death_event <- ifelse(clinical$vital_status == "Alive", 0,1)

# run survival analysis
s <- survfit(Surv(as.numeric(as.character(all_clin$new_death)),all_clin$death_event)~test$status)
s1 <- tryCatch(survdiff(Surv(as.numeric(as.character(all_clin$new_death)),all_clin$death_event)~test$status), error = function(e) return(NA))

# extraect the p.value
pv <- ifelse(is.na(s1),next,(round(1 - pchisq(s1$chisq, length(s1$n) - 1),3)))[[1]]

pv



survdiff(formula = Surv(as.numeric(as.character(all_clin$new_death)),all_clin$death_event)~test$status)

 coxph(Surv(as.numeric(as.character(all_clin$new_death)),all_clin$death_event)~test$status)


# plot the data
plot(survfit(Surv(as.numeric(as.character(all_clin$new_death)),all_clin$death_event)~test$status),
     col=c(1,2), font=2, cex.axis=1.5,lwd=5)

