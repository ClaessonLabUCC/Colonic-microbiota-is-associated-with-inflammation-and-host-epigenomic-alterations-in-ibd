#Author - Feargal Ryan


library(metagenomeSeq)
library(ggplot2)
library(gridExtra)
library(xlsx)
library(gridExtra)

## Reading in tables
mapping = readRDS("SIRG_DADA2_mapping.RDS")
raw_cnts = readRDS("SIRG_DADA2_nc.nh.counts.RDS")
mothur = read.table("SIRG_biopsy_DADA2.mothur.txt",header=TRUE,sep="\t",row.names = 1)
spingo = read.table("SIRG_biopsy_DADA2.final.spingo.7.txt",header=TRUE,sep = "\t",row.names = 1)
family_colors = read.xlsx("family_colors.xlsx",sheetIndex = 1)
rownames(family_colors) = family_colors$family
cnts_filt10 = raw_cnts[apply(raw_cnts>0,1,sum)>=round(ncol(raw_cnts)*0.05),]
cnts_filt10_prop = prop.table(as.matrix(cnts_filt10),2)
cnts_filt10_prop = cnts_filt10_prop * 100

bplots = data.frame(mapping,t(cnts_filt10_prop))


############################
### Random Healthy group ###
############################
## Getting largest sample per patient 
sizes = sort(apply(raw_cnts,2,sum),decreasing=TRUE)
toUnique = data.frame(sizes,"PatientID" = as.character(mapping[names(sizes),"PatientID"]))
toUnique = toUnique[!duplicated(toUnique$PatientID),]
bmap = mapping[rownames(toUnique),]
HTb = rownames(bmap[bmap$Condition=="Healthy",])
mapping$Description3 = mapping$Description2
mapping$Description3[!rownames(mapping) %in% HTb & !mapping$Condition %in% c("CD","UC")] = "HC_n"

##################################
### Create MRexperiment object ###
##################################

## Should there be raw counts filtered
MGS = newMRexperiment(counts = raw_cnts,phenoData = AnnotatedDataFrame(mapping))

## Normalization
MGS_p = cumNormStatFast(MGS)
MGS_css = cumNorm(MGS,p = MGS_p)
norm_cnts = MRcounts(MGS_css,norm=TRUE)
norm_cnts = norm_cnts[rownames(cnts_filt10_prop),]
norm_means = apply(norm_cnts,1,mean)

## Point size calculation
a_s = norm_means
all_size = a_s
all_size[a_s>=as.numeric(summary(a_s)[5])] = "Biggest" # 75th percentile
all_size[a_s>=as.numeric(summary(a_s)[3]) & a_s<as.numeric(summary(a_s)[5]) ] = "Second Biggest"# Median
all_size[a_s>=as.numeric(summary(a_s)[2]) & a_s<as.numeric(summary(a_s)[3]) ] = "Second Smallest"# 25th percentile
all_size[a_s<as.numeric(summary(a_s)[2])] = "Smallest"# Median
all_size = gsub("Second Smallest",5,all_size)
all_size = gsub("Second Biggest",7,all_size)
all_size = gsub("Biggest",10,all_size)
all_size = as.numeric(gsub("Smallest",2,all_size))
names(all_size) = names(a_s)


##################################
###### Differential Testing ######
##################################


## CD inflamed vs Healthy
CDivsHCr = MGS_css[, -which(pData(MGS_css)$Description3 %in% c("CD_NonInflamed","UC_Inflamed","UC_Noninflamed","HC_n") )]
CDivsHCr = filterData(CDivsHCr, present = round(nrow(pData(CDivsHCr))*0.2))
mod <- model.matrix(~ Condition, data = pData(CDivsHCr))
CDivsHCr_fit = fitFeatureModel(obj = CDivsHCr, mod = mod)
CDivsHCr_res = MRfulltable(CDivsHCr_fit,by=2,coef = 2, number = nrow(raw_cnts),eff=.5)
CDivsHCr_sig = CDivsHCr_res[CDivsHCr_res$adjPvalues<0.05,]
CDi_HCr_mothur = data.frame(mothur[rownames(CDivsHCr_res),],spingo[rownames(CDivsHCr_res),"Species"],CDivsHCr_res)
CDi_s = norm_means[rownames(CDivsHCr_res)]
CDivsHCr_res$size = CDi_s
CDivsHCr_res = CDivsHCr_res[order(CDivsHCr_res$size,decreasing=FALSE),]
CDi_sizes = all_size[rownames(CDivsHCr_res)]
CDi_f = as.character(droplevels(mothur[rownames(CDivsHCr_res),"Family"]))
CDi_f = gsub("Desulfovibrionaceae","Other",CDi_f)
CDi_f = gsub("Veillonellaceae","Other",CDi_f)
CDi_f = gsub("Rikenellaceae","Other",CDi_f)
CDi_f = gsub("Peptostreptococcaceae","Other",CDi_f)
CDi_f = gsub("Pasteurellaceae","Other",CDi_f)
CDi_c = as.character(family_colors[levels(as.factor(CDi_f)),"colors"])
CDi_volcano_plot = ggplot(data = CDivsHCr_res,aes(x = CDivsHCr_res$logFC,y = -log10(CDivsHCr_res$adjPvalues),size=as.factor(CDi_sizes),color=CDi_f))+
  geom_point()+geom_hline(yintercept=1.30103)+scale_color_manual("Family",values=CDi_c)+theme_classic()+scale_size_manual(values=c(1,2,5,7))




## CD Noninflamed vs Healthy
CDnivsHCr = MGS_css[, -which(pData(MGS_css)$Description3 %in% c("CD_Inflamed","UC_Inflamed","UC_Noninflamed","HC_n") )]
CDnivsHCr = filterData(CDnivsHCr, present = round(nrow(pData(CDnivsHCr))*0.2),depth=1)
mod <- model.matrix(~ Condition, data = pData(CDnivsHCr))
CDnivsHCr_fit = fitFeatureModel(obj = CDnivsHCr, mod = mod)
CDnivsHCr_res = MRfulltable(CDnivsHCr_fit,by=2,coef = 2, number = nrow(raw_cnts),eff=.5)
CDnivsHCr_sig = CDnivsHCr_res[CDnivsHCr_res$adjPvalues<0.05,]
CDni_HCr_mothur = data.frame(mothur[rownames(CDnivsHCr_res),],spingo[rownames(CDnivsHCr_res),"Species"],CDnivsHCr_res)
CDni_s = norm_means[rownames(CDnivsHCr_res)]
CDnivsHCr_res$size = CDni_s
CDnivsHCr_res = CDnivsHCr_res[order(CDnivsHCr_res$size,decreasing=FALSE),]
CDni_sizes = all_size[rownames(CDnivsHCr_res)]
CDni_sizes[is.na(CDni_sizes)] = 2
CDni_f = as.character(droplevels(mothur[rownames(CDnivsHCr_res),"Family"]))
CDni_f = gsub("Desulfovibrionaceae","Other",CDni_f)
CDni_f = gsub("Veillonellaceae","Other",CDni_f)
CDni_f = gsub("Rikenellaceae","Other",CDni_f)
CDni_f = gsub("Peptostreptococcaceae","Other",CDni_f)
CDni_f = gsub("Pasteurellaceae","Other",CDni_f)
CDni_c = as.character(family_colors[levels(as.factor(CDni_f)),"colors"])
CDni_volcano_plot = ggplot(data = CDnivsHCr_res,aes(x = CDnivsHCr_res$logFC,y = -log10(CDnivsHCr_res$adjPvalues),size=as.factor(CDni_sizes),color=CDni_f))+
  geom_point()+geom_hline(yintercept=1.30103)+scale_color_manual("Family",values=CDni_c)+theme_classic()+scale_size_manual(values=c(1,2,5,7))





## UC inflamed vs Healthy
UCivsHCr = MGS_css[, -which(pData(MGS_css)$Description3 %in% c("UC_Noninflamed","CD_Inflamed","CD_NonInflamed","HC_n") )]
UCivsHCr = filterData(UCivsHCr, present = round(nrow(pData(UCivsHCr))*0.1),depth=1)
mod <- model.matrix(~ Condition, data = pData(UCivsHCr))
UCivsHCr_fit = fitFeatureModel(obj = UCivsHCr, mod = mod)
UCivsHCr_res = MRfulltable(UCivsHCr_fit,by=2,coef = 2, number = nrow(raw_cnts),eff=.5)
UCivsHCr_res = UCivsHCr_res[complete.cases(UCivsHCr_res),]
UCivsHCr_sig = UCivsHCr_res[UCivsHCr_res$adjPvalues<0.05,]
UCi_HCr_mothur = data.frame(mothur[rownames(UCivsHCr_res),],spingo[rownames(UCivsHCr_res),"Species"],UCivsHCr_res)
UCi_s = norm_means[rownames(UCivsHCr_res)]
UCivsHCr_res$size = UCi_s
UCivsHCr_res = UCivsHCr_res[order(UCivsHCr_res$size,decreasing=FALSE),]
UCi_sizes = all_size[rownames(UCivsHCr_res)]
UCi_f = as.character(droplevels(mothur[rownames(UCivsHCr_res),"Family"]))
UCi_f = gsub("Desulfovibrionaceae","Other",UCi_f)
UCi_f = gsub("Veillonellaceae","Other",UCi_f)
UCi_f = gsub("Rikenellaceae","Other",UCi_f)
UCi_f = gsub("Peptostreptococcaceae","Other",UCi_f)
UCi_f = gsub("Pasteurellaceae","Other",UCi_f)
UCi_f = gsub("Acidaminococcaceae","Other",UCi_f)
UCi_f = gsub("Streptococcaceae","Other",UCi_f)
UCi_c = as.character(family_colors[levels(as.factor(UCi_f)),"colors"])
UCi_volcano_plot = ggplot(data = UCivsHCr_res,aes(x = UCivsHCr_res$logFC*-1,y = -log10(UCivsHCr_res$adjPvalues),size=as.factor(UCi_sizes),color=UCi_f))+
  geom_point()+geom_hline(yintercept=1.30103)+scale_color_manual("Family",values=UCi_c)+theme_classic()+scale_size_manual(values=c(1,2,5,7))



## UC Noninflamed vs Healthy
UCnivsHCr = MGS_css[, -which(pData(MGS_css)$Description3 %in% c("UC_Inflamed","CD_Inflamed","CD_NonInflamed","HC_n") )]
UCnivsHCr = filterData(UCnivsHCr, present = round(nrow(pData(UCnivsHCr))*0.1),depth=1)
mod <- model.matrix(~ Condition, data = pData(UCnivsHCr))
UCnivsHCr_fit = fitFeatureModel(obj = UCnivsHCr, mod = mod)
UCnivsHCr_res = MRfulltable(UCnivsHCr_fit,by=2,coef = 2, number = nrow(raw_cnts),eff=.5)
UCnivsHCr_res = UCnivsHCr_res[complete.cases(UCnivsHCr_res),]
UCnivsHCr_sig = UCnivsHCr_res[UCnivsHCr_res$adjPvalues<0.05,]
UCni_HCr_mothur = data.frame(mothur[rownames(UCnivsHCr_res),],spingo[rownames(UCnivsHCr_res),"Species"],UCnivsHCr_res)
UCni_s = norm_means[rownames(UCnivsHCr_res)]
UCnivsHCr_res$size = UCni_s
UCnivsHCr_res = UCnivsHCr_res[order(UCnivsHCr_res$size,decreasing=FALSE),]
UCni_sizes = all_size[rownames(UCnivsHCr_res)]
UCni_f = as.character(droplevels(mothur[rownames(UCnivsHCr_res),"Family"]))
UCni_f = gsub("Desulfovibrionaceae","Other",UCni_f)
UCni_f = gsub("Veillonellaceae","Other",UCni_f)
UCni_f = gsub("Rikenellaceae","Other",UCni_f)
UCni_f = gsub("Peptostreptococcaceae","Other",UCni_f)
UCni_f = gsub("Pasteurellaceae","Other",UCni_f)
UCni_f = gsub("Acidaminococcaceae","Other",UCni_f)
UCni_f = gsub("Streptococcaceae","Other",UCni_f)
UCni_c = as.character(family_colors[levels(as.factor(UCni_f)),"colors"])



UCni_volcano_plot = ggplot(data = UCnivsHCr_res,aes(x = UCnivsHCr_res$logFC*-1,y = -log10(UCnivsHCr_res$adjPvalues),size=as.factor(UCni_sizes),color=UCni_f))+
  geom_point()+geom_hline(yintercept=1.30103)+scale_color_manual("Family",values=UCni_c)+theme_classic()+scale_size_manual(values=c(1,2,5,7))




## CD inflamed vs UC inflamed
CDivsUCi = MGS_css[, -which(pData(MGS_css)$Description3 %in% c("CD_NonInflamed","UC_Noninflamed","HC_n","Healthy") )]
CDivsUCi = filterData(CDivsUCi, present = round(nrow(pData(CDivsUCi))*0.2))
mod <- model.matrix(~ Condition, data = pData(CDivsUCi))
CDivsUCi_fit = fitFeatureModel(obj = CDivsUCi, mod = mod)
CDivsUCi_res = MRfulltable(CDivsUCi_fit,by=2,coef = 2, number = nrow(raw_cnts),eff=.5)
CDivsUCi_res = CDivsUCi_res[complete.cases(CDivsUCi_res),]
CDivsUCi_sig = CDivsUCi_res[CDivsUCi_res$adjPvalues<0.05,]
CDivsUCi_mothur = data.frame(mothur[rownames(CDivsUCi_res),],spingo[rownames(CDivsUCi_res),"Species"],CDivsUCi_res)
CDivsUCi_s = norm_means[rownames(CDivsUCi_res)]
CDivsUCi_res$size = CDivsUCi_s
CDivsUCi_res = CDivsUCi_res[order(CDivsUCi_res$size,decreasing=FALSE),]
CDivsUCi_sizes = all_size[rownames(CDivsUCi_res)]
CDivsUCi_f = as.character(droplevels(mothur[rownames(CDivsUCi_res),"Family"]))
CDivsUCi_f = gsub("Desulfovibrionaceae","Other",CDivsUCi_f)
CDivsUCi_f = gsub("Veillonellaceae","Other",CDivsUCi_f)
CDivsUCi_f = gsub("Rikenellaceae","Other",CDivsUCi_f)
CDivsUCi_f = gsub("Peptostreptococcaceae","Other",CDivsUCi_f)
CDivsUCi_f = gsub("Pasteurellaceae","Other",CDivsUCi_f)
CDivsUCi_f = gsub("Acidaminococcaceae","Other",CDivsUCi_f)
CDivsUCi_f = gsub("Streptococcaceae","Other",CDivsUCi_f)
CDivsUCi_c = as.character(family_colors[levels(as.factor(CDivsUCi_f)),"colors"])
CDivsUCi_volcano_plot = ggplot(data = CDivsUCi_res,aes(x = CDivsUCi_res$logFC,y = -log10(CDivsUCi_res$adjPvalues),size=as.factor(CDivsUCi_sizes),color=CDivsUCi_f))+
  geom_point()+geom_hline(yintercept=1.30103)+scale_color_manual("Family",values=CDivsUCi_c)+theme_classic()+scale_size_manual(values=c(1,2,5,7))


## CD noninflamed vs UC noninflamed
CDnivsUCni = MGS_css[, -which(pData(MGS_css)$Description3 %in% c("CD_Inflamed","UC_Inflamed","HC_n","Healthy") )]
CDnivsUCni = filterData(CDnivsUCni, present = round(nrow(pData(CDnivsUCni))*0.2))
mod <- model.matrix(~ Condition, data = pData(CDnivsUCni))
CDnivsUCni_fit = fitFeatureModel(obj = CDnivsUCni, mod = mod)
CDnivsUCni_res = MRfulltable(CDnivsUCni_fit,by=2,coef = 2, number = nrow(raw_cnts),eff=.5)
CDnivsUCni_res = CDnivsUCni_res[complete.cases(CDnivsUCni_res),]
CDnivsUCni_sig = CDnivsUCni_res[CDnivsUCni_res$adjPvalues<0.05,]
CDnivsUCni_mothur = data.frame(mothur[rownames(CDnivsUCni_res),],spingo[rownames(CDnivsUCni_res),"Species"],CDnivsUCni_res)
CDnivsUCni_s = norm_means[rownames(CDnivsUCni_res)]
CDnivsUCni_res$size = CDnivsUCni_s
CDnivsUCni_res = CDnivsUCni_res[order(CDnivsUCni_res$size,decreasing=FALSE),]
CDnivsUCni_sizes = all_size[rownames(CDnivsUCni_res)]
CDnivsUCni_f = as.character(droplevels(mothur[rownames(CDnivsUCni_res),"Family"]))
CDnivsUCni_f = gsub("Desulfovibrionaceae","Other",CDnivsUCni_f)
CDnivsUCni_f = gsub("Veillonellaceae","Other",CDnivsUCni_f)
CDnivsUCni_f = gsub("Rikenellaceae","Other",CDnivsUCni_f)
CDnivsUCni_f = gsub("Peptostreptococcaceae","Other",CDnivsUCni_f)
CDnivsUCni_f = gsub("Pasteurellaceae","Other",CDnivsUCni_f)
CDnivsUCni_f = gsub("Acidaminococcaceae","Other",CDnivsUCni_f)
CDnivsUCni_f = gsub("Streptococcaceae","Other",CDnivsUCni_f)
table(levels(as.factor(CDnivsUCni_f)) %in% family_colors$family)
CDnivsUCni_c = as.character(family_colors[levels(as.factor(CDnivsUCni_f)),"colors"])
CDnivsUCni_volcano_plot = ggplot(data = CDnivsUCni_res,aes(x = CDnivsUCni_res$logFC,y = -log10(CDnivsUCni_res$adjPvalues),size=as.factor(CDnivsUCni_sizes),color=CDnivsUCni_f))+
  geom_point()+geom_hline(yintercept=1.30103)+scale_color_manual("Family",values=CDnivsUCni_c)+theme_classic()+scale_size_manual(values=c(1,2,5,7))








###############################################################################################################
###############################################################################################################

paired_mapping = mapping[mapping$Paired=="Yes",]
paired_cnts = raw_cnts[,rownames(paired_mapping)]
pMGS = newMRexperiment(counts = paired_cnts,phenoData = AnnotatedDataFrame(paired_mapping))

## Normalization
pMGS_p = cumNormStatFast(pMGS)
pMGS_css = cumNorm(pMGS,p = pMGS_p)
pnorm_cnts = MRcounts(pMGS_css,norm=TRUE)
pnorm_means = apply(pnorm_cnts,1,mean)


## CD inflamed vs CD noninflamed
CDnivsCDi = pMGS_css[, -which(pData(pMGS_css)$Description3 %in% c("UC_Noninflamed","UC_Inflamed","HC_n","Healthy") )]
CDnivsCDi = filterData(CDnivsCDi, present = round(nrow(pData(CDnivsCDi))*0.2))
mod <- model.matrix(~ Status, data = pData(CDnivsCDi))
CDnivsCDi_fit = fitFeatureModel(obj = CDnivsCDi, mod = mod)
CDnivsCDi_res = MRfulltable(CDnivsCDi_fit,by=2,coef = 2, number = nrow(raw_cnts),eff=.5)
CDnivsCDi_res = CDnivsCDi_res[complete.cases(CDnivsCDi_res),]
CDnivsCDi_mothur = data.frame(mothur[rownames(CDnivsCDi_res),],spingo[rownames(CDnivsCDi_res),"Species"],CDnivsCDi_res)
CDnivsCDi_sig = CDnivsCDi_res[CDnivsCDi_res$adjPvalues<0.05,]
CDnivsCDi_s = norm_means[rownames(CDnivsCDi_res)]
CDnivsCDi_res$size = CDnivsCDi_s
CDnivsCDi_res = CDnivsCDi_res[order(CDnivsCDi_res$size,decreasing=FALSE),]
CDnivsCDi_sizes = all_size[rownames(CDnivsCDi_res)]
CDnivsCDi_f = as.character(droplevels(mothur[rownames(CDnivsCDi_res),"Family"]))
CDnivsCDi_f = gsub("Desulfovibrionaceae","Other",CDnivsCDi_f)
CDnivsCDi_f = gsub("Veillonellaceae","Other",CDnivsCDi_f)
CDnivsCDi_f = gsub("Rikenellaceae","Other",CDnivsCDi_f)
CDnivsCDi_f = gsub("Peptostreptococcaceae","Other",CDnivsCDi_f)
CDnivsCDi_f = gsub("Pasteurellaceae","Other",CDnivsCDi_f)
CDnivsCDi_f = gsub("Acidaminococcaceae","Other",CDnivsCDi_f)
CDnivsCDi_f = gsub("Streptococcaceae","Other",CDnivsCDi_f)
CDnivsCDi_c = as.character(family_colors[levels(as.factor(CDnivsCDi_f)),"colors"])
CDnivsCDi_volcano_plot = ggplot(data = CDnivsCDi_res,aes(x = CDnivsCDi_res$logFC,y = -log10(CDnivsCDi_res$adjPvalues),size=as.factor(CDnivsCDi_sizes),color=CDnivsCDi_f))+
  geom_point()+geom_hline(yintercept=1.30103)+scale_color_manual("Family",values=CDnivsCDi_c)+theme_classic()+scale_size_manual(values=c(1,2,5,7))







## UC inflamed vs UC noninflamed
UCnivsUCi = pMGS_css[, -which(pData(pMGS_css)$Condition %in% c("CD","Healthy") )]
UCnivsUCi = filterData(UCnivsUCi, present = round(nrow(pData(UCnivsUCi))*0.2))
mod <- model.matrix(~ Status, data = pData(UCnivsUCi))
UCnivsUCi_fit = fitFeatureModel(obj = UCnivsUCi, mod = mod)
UCnivsUCi_res = MRfulltable(UCnivsUCi_fit,by=2,coef = 2, number = nrow(raw_cnts),eff=.5)
UCnivsUCi_res = UCnivsUCi_res[complete.cases(UCnivsUCi_res),]
UCnivsUCi_mothur = data.frame(mothur[rownames(UCnivsUCi_res),],spingo[rownames(UCnivsUCi_res),"Species"],UCnivsUCi_res)
UCnivsUCi_sig = UCnivsUCi_res[UCnivsUCi_res$adjPvalues<0.05,]
UCnivsUCi_s = norm_means[rownames(UCnivsUCi_res)]
UCnivsUCi_res$size = UCnivsUCi_s
UCnivsUCi_res = UCnivsUCi_res[order(UCnivsUCi_res$size,decreasing=FALSE),]
UCnivsUCi_sizes = all_size[rownames(UCnivsUCi_res)]
UCnivsUCi_f = as.character(droplevels(mothur[rownames(UCnivsUCi_res),"Family"]))
UCnivsUCi_f = gsub("Desulfovibrionaceae","Other",UCnivsUCi_f)
UCnivsUCi_f = gsub("Veillonellaceae","Other",UCnivsUCi_f)
UCnivsUCi_f = gsub("Rikenellaceae","Other",UCnivsUCi_f)
UCnivsUCi_f = gsub("Peptostreptococcaceae","Other",UCnivsUCi_f)
UCnivsUCi_f = gsub("Pasteurellaceae","Other",UCnivsUCi_f)
UCnivsUCi_f = gsub("Acidaminococcaceae","Other",UCnivsUCi_f)
UCnivsUCi_f = gsub("Streptococcaceae","Other",UCnivsUCi_f)
UCnivsUCi_c = as.character(family_colors[levels(as.factor(UCnivsUCi_f)),"colors"])
UCnivsUCi_volcano_plot = ggplot(data = UCnivsUCi_res,aes(x = UCnivsUCi_res$logFC,y = -log10(UCnivsUCi_res$adjPvalues),size=as.factor(UCnivsUCi_sizes),color=UCnivsUCi_f))+
  geom_point()+geom_hline(yintercept=1.30103)+scale_color_manual("Family",values=UCnivsUCi_c)+theme_classic()+scale_size_manual(values=c(1,2,5,7))





#### Combined plot of all volcano plots ####
grid.arrange(
  CDnivsCDi_volcano_plot+theme(axis.title.x=element_blank(),legend.position="none",axis.title.y=element_blank())+scale_y_continuous(limits=c(0,4.5))+scale_x_continuous(limits=c(-2.5,2.5)),
  CDi_volcano_plot+theme(axis.title.x=element_blank(),legend.position="none",axis.title.y=element_blank())+scale_y_continuous(limits=c(0,4.5))+scale_x_continuous(limits=c(-2.5,2.5)),
  CDni_volcano_plot+theme(axis.title.x=element_blank(),legend.position="none",axis.title.y=element_blank())+scale_y_continuous(limits=c(0,4.5))+scale_x_continuous(limits=c(-2.5,2.5)),
  CDivsUCi_volcano_plot+theme(axis.title.x=element_blank(),legend.position="none",axis.title.y=element_blank())+scale_y_continuous(limits=c(0,4.5))+scale_x_continuous(limits=c(-2.5,2.5)),
  UCnivsUCi_volcano_plot+theme(axis.title.x=element_blank(),legend.position="none",axis.title.y=element_blank())+scale_y_continuous(limits=c(0,4.5))+scale_x_continuous(limits=c(-2.5,2.5)),
  UCi_volcano_plot+theme(axis.title.x=element_blank(),legend.position="none",axis.title.y=element_blank())+scale_y_continuous(limits=c(0,4.5))+scale_x_continuous(limits=c(-2.5,2.5)),
  UCni_volcano_plot+theme(axis.title.x=element_blank(),legend.position="none",axis.title.y=element_blank())+scale_y_continuous(limits=c(0,4.5))+scale_x_continuous(limits=c(-2.5,2.5)),
  CDnivsUCni_volcano_plot+theme(axis.title.x=element_blank(),legend.position="none",axis.title.y=element_blank())+scale_y_continuous(limits=c(0,4.5))+scale_x_continuous(limits=c(-2.5,2.5)),
nrow=2)
