if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("RPMM")



library("minfi")
#R code to convert IDAT files to raw beta values
dir="~/PhD/Data/Methylation/meth"  #working directory
#must be only file in folder

target=read.metharray.sheet(dir) #Read in samplesheet in csv format
RGset=read.metharray.exp(targets=target)

pd=pData(RGset)  #Access Pheno date
#target=read.metharray.sheet(dir) #Read in samplesheet in csv format
#dir="~/PhD/Data/Methylation/"  #working directory
#target=read.metharray.sheet(dir) #Read in samplesheet in csv format
#dir="~/PhD/Data/Methylation"  #working directory
#target=read.metharray.sheet(dir) #Read in samplesheet in csv format
#RGset=read.metharray.exp(targets=target)
#pd=pData(RGset)  #Access Pheno date
#Preprocessing
#raw and betas
MSet.raw=preprocessRaw(RGset)
bvalraw=getBeta(MSet.raw)
#Failed probes
pdet=detectionP(RGset,type="m+u")
failed=pdet>0.01
failedProbes=rownames(failed)[rowMeans(failed)>0.01]
#Filter out failedProbes
raw=bvalraw[!rownames(bvalraw)%in%failedProbes,]

XYchrome <- grep("x|y", rownames(raw), ignore.case = TRUE)
rawT <- raw[XYchrome,]
rows.to.lose<-rownames(rawT)
'%!in%' <- function(x,y)!('%in%'(x,y))

raw=raw[!rownames(raw) %in% rows.to.lose,]
############################################################
mani=read.csv("~/PhD/HumanMethylation450_15017482_v1-2.csv",header=TRUE)
lumi=mani[,c(1,7,12,13,22,23)]

#remove SNPs from data
snplist=lumi[grep('rs',lumi[,1]),]
#Remove probes that multiple align to various sites
probe_file=read.table("~/PhD/Data/Methylation/Probe_info_17764_bad_header.txt",header=TRUE)
probe=probe_file$cg

raw=raw[!rownames(raw)%in%probe,]

#Annotate with type and cg identifiers
lumi = mani[mani$IlmnID%in%rownames(raw),c("IlmnID","Infinium_Design_Type")]
#raw1 <- as.data.frame(raw)
#raw1$IlmnID <- rownames(raw1)
#library(plyr)
#rawie = join(raw1,lumi, by="IlmnID") 
rawie = merge(lumi,raw,by.x="IlmnID",by.y="row.names")  #merge lumi&raw


#Density Plot
library(preprocessCore)

rowN <- rawie[,1]
labels=rawie[,2]
labels=gsub("II","Type 2",labels)
labels=gsub("I","Type 1",labels)

densityPlot(as.matrix(rawie[,-c(1:2)]),sampGroups=labels,main="Beta Values",pal=c("red","plum2"))

#Needed for BMIQ
library(RPMM)
library(methylumi)
library(wateRmelon)

# BMIQ  for your raw betas and then qqnorm. Then prcomp/princomp for PCA. 
RawA <- rawie

RawA$Infinium_Design_Type <- as.character(RawA$Infinium_Design_Type)


RawA$Infinium_Design_Type <- replace(RawA$Infinium_Design_Type , RawA$Infinium_Design_Type == "I", 1)

RawA$Infinium_Design_Type <- factor(RawA$Infinium_Design_Type)
table(RawA$Infinium_Design_Type)
nlevels(RawA$Infinium_Design_Type)
RawA$Infinium_Design_Type <- as.integer(RawA$Infinium_Design_Type)

#normBeta <- data.frame(matrix(, nrow=483970, ncol=108))

normBeta <- data.frame(matrix(, nrow=466379, ncol=108))
for (i in (1:108)){
  j = i+2
  test2 <- BMIQ(RawA[,j], RawA$Infinium_Design_Type)
  normBeta[,i] <-test2$nbeta
}

colnames(normBeta) <- colnames(RawA[3:110])
normBeta$IlmnID <- RawA$IlmnID



# normalising the probes to N(0,1)
normalize <- function(x) {
  n <- qqnorm(x, plot.it=FALSE)
  return(n$x)
  
}

#qq <- (apply(normBeta[,-c(109)],1,normalize))
qq <- t(apply(normBeta[,-c(109)],1,normalize)) #Check that you need to transpose here
# qq1 <- (apply(normBeta[,-c(109)],1,normalize))
#Here is the code to check NA's and use only complete cases. 
anyNA(qq)
colSums(is.na(qq))
colnames(qq) <- colnames(normBeta[,-c(109)])
rownames(qq) <- normBeta$IlmnID
qq.ncc=qq[!complete.cases(qq),]#Incomplete cases, which have missing data
qq.ready=qq[complete.cases(qq),]#complete cases


##############
lumi.ready <- lumi[complete.cases(qq),]

metaD <- read.csv("~/PhD/Data/Methylation/meth/metharray_Samplesheet.csv")
metaD$SampleID.1 <- metaD$Sample_Name

colnames(qq) <- metaD$Sample_Name
metaD$SampleID.1
for (i in 1:108){print(colnames(qq.ready)[i]) 
  print(paste(metaD$Sentrix_ID[i], metaD$Sentrix_Position[i], sep = "_"))}

qq.readyBKUP <- qq.ready

colnames(qq.readyBKUP) <- metaD$SampleID.1
colnames(qq.ready)
for (i in 1:108){print(colnames(qq.ready)[i])
  print(paste(metaD$Sentrix_ID[i], metaD$Sentrix_Position[i], sep = "_"))}
qq.readyBKUP <- qq.ready
colnames(qq.readyBKUP) <- metaD$SampleID.1
for (i in 1:108){print(colnames(qq.ready)[i])
  print(paste(metaD$Sentrix_ID[i], metaD$Sentrix_Position[i], sep = "_"))}
colnames(qq.ready) <- metaD$SampleID.1
colnames(qq.ready)
no16s <- c("sCD28i"  , "sCD55a"  , "sCD50a_B", "sCD15a" ,  "sCD15i"   ,"sCD44a" ,  "sCD02a"
)
colQQR <- colnames(qq.ready)

intersect(colQQR, no16s)
setdiff(colQQR, no16s)
keep <- setdiff(colQQR, no16s)
keep[76]
keep <- keep[-76]

qq.ready <- qq.ready[,keep]

#The PCA ustilises the CpG sites values only, so you can keep the normalised beta values separate from the meta data.
#PCA using prcomp, also look into princomp as this can also be used depending on your data
pcareresults=prcomp(t(qq.ready))

#Plot to show variance for first  10 components
results=round(((pcareresults$sdev)^2 / sum(pcareresults$sdev^2))*100,digits=2)
labels=as.character(results[1:10])
labels=paste(labels,"%",sep="")
barplot(results[1:10], names=labels,ylim=c(0,25 ),col="yellow")


#ggplot2 for PC1 vs PC2, colouring for condition, shapes for status and lines to connect samples from same patient.
labx1=paste("PC1",labels[1],sep=" ")
laby1=paste("PC2",labels[2],sep=" ")
laby2=paste("PC4",labels[4],sep=" ")
lab3=paste("PC3",labels[3],sep=" ")
lab5=paste("PC5",labels[5],sep=" ")
lab6=paste("PC6",labels[6],sep=" ")
lab7=paste("PC7",labels[7],sep=" ")
lab8=paste("PC8",labels[8],sep=" ")
lab9=paste("PC9",labels[9],sep=" ")
lab10=paste("PC10",labels[10],sep=" ")
pc1=pcareresults$x[,1]
pc2=pcareresults$x[,2]

pc3=pcareresults$x[,3]
pc4=pcareresults$x[,4]

pc5=pcareresults$x[,5]
pc6=pcareresults$x[,6]

pc7=pcareresults$x[,7]
pc8=pcareresults$x[,8]

pc9=pcareresults$x[,9]
pc10=pcareresults$x[,10]

