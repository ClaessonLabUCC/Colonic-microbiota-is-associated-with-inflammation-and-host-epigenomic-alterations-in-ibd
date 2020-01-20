# raw betas boxplots

library(ggplot2)
library(ggbeeswarm)

library("minfi")
#R code to convert IDAT files to raw beta values
dir="~/PhD/Data/Methylation/meth"  #working directory
#must be only file in folder

target=read.metharray.sheet(dir) #Read in samplesheet in csv format
RGset=read.metharray.exp(targets=target)

pd=pData(RGset)  #Access Pheno date

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

RawA <- rawie

RawA$Infinium_Design_Type <- as.character(RawA$Infinium_Design_Type)


RawA$Infinium_Design_Type <- replace(RawA$Infinium_Design_Type , RawA$Infinium_Design_Type == "I", 1)

RawA$Infinium_Design_Type <- factor(RawA$Infinium_Design_Type)
table(RawA$Infinium_Design_Type)
nlevels(RawA$Infinium_Design_Type)
RawA$Infinium_Design_Type <- as.integer(RawA$Infinium_Design_Type)
rownames(RawA) <- RawA$IlmnID

rawBetas <- RawA[3:110]
dim(rawBetas)


metaD <- read.csv("metharray_Samplesheet.csv")

for (i in 1:108){print(colnames(rawBetas)[i]) 
  print(paste(metaD$Sentrix_ID[i], metaD$Sentrix_Position[i], sep = "_"))}

colnames(rawBetas) <- metaD$Sample_Name


load("100_samples_updated_metadata.RData")


rawBeta100 <- rawBetas[,meta100$SampleID.1]
dim(rawBeta100)
length(intersect(colnames(rawBeta100), meta100$SampleID.1))


for (i in 1:100){print(colnames(rawBeta100)[i]) 
  print(meta100$SampleID.1[i])}

#cg00000321
head(rownames(rawBeta100))

rawBeta100Flip <- t(rawBeta100)


##############################
######### CLUSTERS ###########
##############################

### CLUSTER 1 ###
#cg20044211 NOTCH4

RBFlipcg20044211<- rawBeta100Flip[,c("cg20044211")]

meta100$cg20044211 <- RBFlipcg20044211 

meta100temp <- meta100

metaAll_temp <- meta100temp
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)
metaAll_temp$cluster[metaAll_temp$cluster != "1"] <- "Other"
metaAll_temp$cluster <- as.factor(metaAll_temp$cluster)
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)

ggplot(metaAll_temp, aes(x=cluster, y=cg20044211, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("NOTCH4 (cg20044211)") +  scale_fill_manual("cluster",values=c("#CC0066","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))

### CLUSTER2 ###
# cg07884764 CCDC88B

RBFlipcg07884764<- rawBeta100Flip[,c("cg07884764")]

meta100$cg07884764 <- RBFlipcg07884764 

meta100temp <- meta100

metaAll_temp <- meta100temp
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)
metaAll_temp$cluster[metaAll_temp$cluster != "2"] <- "Other"
metaAll_temp$cluster <- as.factor(metaAll_temp$cluster)
class(metaAll_temp$cluster)

ggplot(metaAll_temp, aes(x=cluster, y=cg07884764, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("CCDC88B (cg07884764)") +  scale_fill_manual("cluster",values=c("#27F0C9","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))


# cg19284131 TRIM27

RBFlipcg19284131<- rawBeta100Flip[,c("cg19284131")]

meta100$cg19284131 <- RBFlipcg19284131 

meta100temp <- meta100

metaAll_temp <- meta100temp
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)
metaAll_temp$cluster[metaAll_temp$cluster != "2"] <- "Other"
metaAll_temp$cluster <- as.factor(metaAll_temp$cluster)
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)

ggplot(metaAll_temp, aes(x=cluster, y=cg19284131, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("TRIM27 (cg19284131)") +  scale_fill_manual("cluster",values=c("#27C0C1","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))


# cg18045172 TAP2

RBFlipcg18045172<- rawBeta100Flip[,c("cg18045172")]

meta100$cg18045172 <- RBFlipcg18045172 

meta100temp <- meta100

metaAll_temp <- meta100temp
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)
metaAll_temp$cluster[metaAll_temp$cluster != "2"] <- "Other"
metaAll_temp$cluster <- as.factor(metaAll_temp$cluster)
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)

ggplot(metaAll_temp, aes(x=cluster, y=cg18045172, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("TAP2 (cg18045172)") +  scale_fill_manual("cluster",values=c("#FFC0F1","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))

### CLUSTER 3 ###


# 	cg02754643 DRAM1
RBFlipcg02754643<- rawBeta100Flip[,c("cg02754643")]

meta100$cg02754643 <- RBFlipcg02754643 

meta100temp <- meta100

metaAll_temp <- meta100temp
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)
metaAll_temp$cluster[metaAll_temp$cluster != "3"] <- "Other"
metaAll_temp$cluster <- as.factor(metaAll_temp$cluster)
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)

ggplot(metaAll_temp, aes(x=cluster, y=cg02754643, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("DRAM1 (cg02754643)") +  scale_fill_manual("cluster",values=c("#F47B73","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))



#########################################
### inflamed vs non inflamed
########################################


# cg03211132 SELE

metaCD <- meta100[meta100$Condition =="CD",]

rawBeta100FlipCD <- rawBeta100Flip[metaCD$SampleID.1,]

RBFlipcg03211132<- rawBeta100FlipCD[,c("cg03211132")]

metaCD$cg03211132 <- RBFlipcg03211132

meta100temp <- metaCD

metaAll_temp <- meta100temp
class(metaAll_temp$Status)
table(metaAll_temp$Status)

ggplot(metaAll_temp, aes(x=Status, y=cg03211132, fill=Status)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("SELE (cg03211132)") +  scale_fill_manual("cluster",values=c("#FFFFFF","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))



##########################################
### SUPPLEMENTARY FIGURE CLUSTER PLOTS ###
##########################################

## CLUSTER 1

# PODXL2	cg08658371

RBFlipcg08658371<- rawBeta100Flip[,c("cg08658371")]

meta100$cg08658371 <- RBFlipcg08658371 

meta100temp <- meta100

metaAll_temp <- meta100temp
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)
metaAll_temp$cluster[metaAll_temp$cluster != "1"] <- "Other"
metaAll_temp$cluster <- as.factor(metaAll_temp$cluster)
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)

ggplot(metaAll_temp, aes(x=cluster, y=cg08658371, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("PODXL2	 (cg08658371)") +  scale_fill_manual("cluster",values=c("#CC0066","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))

C1plot1 <-ggplot(metaAll_temp, aes(x=cluster, y=cg08658371, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("PODXL2	 (cg08658371)") +  scale_fill_manual("cluster",values=c("#CC0066","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))

# OPN5	cg11801512

RBFlipcg11801512<- rawBeta100Flip[,c("cg11801512")]

meta100$cg11801512 <- RBFlipcg11801512 

meta100temp <- meta100

metaAll_temp <- meta100temp
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)
metaAll_temp$cluster[metaAll_temp$cluster != "1"] <- "Other"
metaAll_temp$cluster <- as.factor(metaAll_temp$cluster)
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)

ggplot(metaAll_temp, aes(x=cluster, y=cg11801512, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("OPN5	 (cg11801512)") +  scale_fill_manual("cluster",values=c("#CC0066","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))

C1plot2 <- ggplot(metaAll_temp, aes(x=cluster, y=cg11801512, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("OPN5	 (cg11801512)") +  scale_fill_manual("cluster",values=c("#CC0066","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))



# cg01732952 RAVER2

RBFlipcg01732952<- rawBeta100Flip[,c("cg01732952")]

meta100$cg01732952 <- RBFlipcg01732952 

meta100temp <- meta100

metaAll_temp <- meta100temp
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)
metaAll_temp$cluster[metaAll_temp$cluster != "1"] <- "Other"
metaAll_temp$cluster <- as.factor(metaAll_temp$cluster)
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)

ggplot(metaAll_temp, aes(x=cluster, y= cg01732952, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("RAVER2	 (cg01732952)") +  scale_fill_manual("cluster",values=c("#CC0066","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))

C1plot3 <- ggplot(metaAll_temp, aes(x=cluster, y= cg01732952, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("RAVER2	 (cg01732952)") +  scale_fill_manual("cluster",values=c("#CC0066","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))


## CLUSTER 2

#  LOXL2	cg09042448

RBFlipcg09042448<- rawBeta100Flip[,c("cg09042448")]

meta100$cg09042448 <- RBFlipcg09042448 

meta100temp <- meta100

metaAll_temp <- meta100temp
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)
metaAll_temp$cluster[metaAll_temp$cluster != "2"] <- "Other"
metaAll_temp$cluster <- as.factor(metaAll_temp$cluster)
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)

ggplot(metaAll_temp, aes(x=cluster, y= cg09042448, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("LOXL2	 (cg09042448)") +  scale_fill_manual("cluster",values=c("#663399","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))

C2plot1 <- ggplot(metaAll_temp, aes(x=cluster, y= cg09042448, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("LOXL2	 (cg09042448)") +  scale_fill_manual("cluster",values=c("#663399","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))

# RAB33A	cg05911704

RBFlipcg05911704<- rawBeta100Flip[,c("cg05911704")]

meta100$cg05911704 <- RBFlipcg05911704 

meta100temp <- meta100

metaAll_temp <- meta100temp
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)
metaAll_temp$cluster[metaAll_temp$cluster != "2"] <- "Other"
metaAll_temp$cluster <- as.factor(metaAll_temp$cluster)
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)

ggplot(metaAll_temp, aes(x=cluster, y= cg05911704, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("RAB33A	 (cg05911704)") +  scale_fill_manual("cluster",values=c("#663399","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))

C2plot2 <- ggplot(metaAll_temp, aes(x=cluster, y= cg05911704, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("RAB33A	 (cg05911704)") +  scale_fill_manual("cluster",values=c("#663399","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))

ggplot(metaAll_temp, aes(x=cluster, y= cg05911704, fill=cluster,colour=Status)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("RAB33A	 (cg05911704)") +  scale_fill_manual("cluster",values=c("#663399","#808080")) + 
  theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))

ggplot(metaAll_temp, aes(x=cluster, y= cg05911704, fill=cluster,colour=Sex)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("RAB33A	 (cg05911704)") +  scale_fill_manual("cluster",values=c("#663399","#808080")) + 
  theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))

metaAll_temp$Sex <- as.character(metaAll_temp$Sex)
XmetaAll_temp <- metaAll_temp[metaAll_temp$Sex == "Female",]
C2plot2 <- ggplot(XmetaAll_temp, aes(x=cluster, y= cg05911704, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("RAB33A	 (cg05911704)") +  scale_fill_manual("cluster",values=c("#663399","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))

ggplot(metaAll_temp, aes(x=cluster, y= cg05911704, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3,size=1)+
  ggtitle("RAB33A	 (cg05911704)") +  scale_fill_manual("cluster",values=c("#663399","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))


# ATP4A	cg02456521

RBFlipcg02456521<- rawBeta100Flip[,c("cg02456521")]

meta100$cg02456521 <- RBFlipcg02456521 

meta100temp <- meta100

metaAll_temp <- meta100temp
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)
metaAll_temp$cluster[metaAll_temp$cluster != "2"] <- "Other"
metaAll_temp$cluster <- as.factor(metaAll_temp$cluster)
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)

ggplot(metaAll_temp, aes(x=cluster, y= cg02456521, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("ATP4A	 (cg02456521)") +  scale_fill_manual("cluster",values=c("#663399","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))

C2plot3 <- ggplot(metaAll_temp, aes(x=cluster, y= cg02456521, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("ATP4A	 (cg02456521)") +  scale_fill_manual("cluster",values=c("#663399","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))


## CLUSTER 3

# KCNQ1	cg19728223

RBFlipcg19728223<- rawBeta100Flip[,c("cg19728223")]

meta100$cg19728223 <- RBFlipcg19728223 

meta100temp <- meta100

metaAll_temp <- meta100temp
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)
metaAll_temp$cluster[metaAll_temp$cluster != "3"] <- "Other"
metaAll_temp$cluster <- as.factor(metaAll_temp$cluster)
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)

ggplot(metaAll_temp, aes(x=cluster, y= cg19728223, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("KCNQ1	 (cg19728223)") +  scale_fill_manual("cluster",values=c("#6666FF","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))

C3plot1 <- ggplot(metaAll_temp, aes(x=cluster, y= cg19728223, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("KCNQ1	 (cg19728223)") +  scale_fill_manual("cluster",values=c("#6666FF","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))


# SOX2OT	cg25960893

RBFlipcg25960893<- rawBeta100Flip[,c("cg25960893")]

meta100$cg25960893 <- RBFlipcg25960893 

meta100temp <- meta100

metaAll_temp <- meta100temp
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)
metaAll_temp$cluster[metaAll_temp$cluster != "3"] <- "Other"
metaAll_temp$cluster <- as.factor(metaAll_temp$cluster)
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)

ggplot(metaAll_temp, aes(x=cluster, y= cg25960893, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("SOX2OT	 (cg25960893)") +  scale_fill_manual("cluster",values=c("#6666FF","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))

C3plot2 <- ggplot(metaAll_temp, aes(x=cluster, y= cg25960893, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("SOX2OT	 (cg25960893)") +  scale_fill_manual("cluster",values=c("#6666FF","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))


# LHX5	cg23497704

RBFlipcg23497704<- rawBeta100Flip[,c("cg23497704")]

meta100$cg23497704 <- RBFlipcg23497704 

meta100temp <- meta100

metaAll_temp <- meta100temp
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)
metaAll_temp$cluster[metaAll_temp$cluster != "3"] <- "Other"
metaAll_temp$cluster <- as.factor(metaAll_temp$cluster)
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)

ggplot(metaAll_temp, aes(x=cluster, y= cg23497704, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("LHX5	 (cg23497704)") +  scale_fill_manual("cluster",values=c("#6666FF","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))

C3plot3 <- ggplot(metaAll_temp, aes(x=cluster, y= cg23497704, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("LHX5	 (cg23497704)") +  scale_fill_manual("cluster",values=c("#6666FF","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))


## CLUSTER 4

# TNK2	cg27418204


RBFlipcg27418204<- rawBeta100Flip[,c("cg27418204")]

meta100$cg27418204 <- RBFlipcg27418204 

meta100temp <- meta100

metaAll_temp <- meta100temp
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)
metaAll_temp$cluster[metaAll_temp$cluster != "4"] <- "Other"
metaAll_temp$cluster <- as.factor(metaAll_temp$cluster)
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)

ggplot(metaAll_temp, aes(x=cluster, y= cg27418204, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("TNK2	 (cg27418204)") +  scale_fill_manual("cluster",values=c("#66CCFF","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))

C4plot1 <- ggplot(metaAll_temp, aes(x=cluster, y= cg27418204, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("TNK2	 (cg27418204)") +  scale_fill_manual("cluster",values=c("#66CCFF","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))


# ESD	cg11946072

RBFlipcg11946072<- rawBeta100Flip[,c("cg11946072")]

meta100$cg11946072 <- RBFlipcg11946072 

meta100temp <- meta100

metaAll_temp <- meta100temp
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)
metaAll_temp$cluster[metaAll_temp$cluster != "4"] <- "Other"
metaAll_temp$cluster <- as.factor(metaAll_temp$cluster)
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)

ggplot(metaAll_temp, aes(x=cluster, y= cg11946072, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("ESD	 (cg11946072)") +  scale_fill_manual("cluster",values=c("#66CCFF","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))

C4plot2 <- ggplot(metaAll_temp, aes(x=cluster, y= cg11946072, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("ESD	 (cg11946072)") +  scale_fill_manual("cluster",values=c("#66CCFF","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))


# ARIH2	cg14039083
RBFlipcg14039083<- rawBeta100Flip[,c("cg14039083")]

meta100$cg14039083 <- RBFlipcg14039083 

meta100temp <- meta100

metaAll_temp <- meta100temp
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)
metaAll_temp$cluster[metaAll_temp$cluster != "4"] <- "Other"
metaAll_temp$cluster <- as.factor(metaAll_temp$cluster)
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)

ggplot(metaAll_temp, aes(x=cluster, y= cg14039083, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("ARIH2	 (cg14039083)") +  scale_fill_manual("cluster",values=c("#66CCFF","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))

C4plot3 <- ggplot(metaAll_temp, aes(x=cluster, y= cg14039083, fill=cluster)) +geom_boxplot(outlier.shape = NA) + geom_beeswarm(cex = 3)+
  ggtitle("ARIH2	 (cg14039083)") +  scale_fill_manual("cluster",values=c("#66CCFF","#808080")) + 
  theme(legend.position = "none",panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA), plot.title = element_text(hjust = 0.5))

library(gridExtra)

grid.arrange(C1plot1,C1plot2, C1plot3, C2plot1,C2plot2, C2plot3,C3plot1,C3plot2, C3plot3,C4plot1,C4plot2, C4plot3)
