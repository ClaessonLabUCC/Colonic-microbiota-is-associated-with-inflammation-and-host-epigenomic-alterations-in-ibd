# load in normalised betas -> qq.ready
# loading in metadata for 100 CD samples -> meta100


load("100_samples_updated_metadata.RData")
library(lme4)

LinMod466209 <- data.frame(matrix(, nrow=466209, ncol=2))

metaAll_temp <- meta100
qqDF_temp1 <- qq.ready

metaAll_temp <- metaAll_temp[metaAll_temp$Condition == "CD",]
qqDF_temp1 <- qqDF_temp1[,metaAll_temp$SampleID.1]
  
class(metaAll_temp$Status)
table(metaAll_temp$Status)


# Run linear mixed model

for (i in (1:466209)){
  lmmT <- NULL
  lmm2T <- NULL
  print(rownames(qqDF_temp1)[i])
  x <- rownames(qqDF_temp1)[i]
  metaAll_temp1 <- metaAll_temp
  metaAll_temp1$temp <- t(qqDF_temp1[i,])[1,]
  #print(t(qqDF_temp1[i,]))
  lmmT <- try(lmer(temp ~ Status + Sex +  Smoker + Age + Sentrix_ID + Sentrix_Position  + Biopsy_location +Hospital+ (1 | Patient_ID), data = metaAll_temp1,
                   REML = FALSE))
  lmm2T <- try(lmer(temp ~ Sex + Smoker + Age + Sentrix_ID + Sentrix_Position  + Biopsy_location + Hospital+(1 | Patient_ID), data = metaAll_temp1,
                    REML = FALSE))
  #print(summary(lmm))
  #print(anova(lmm2T, lmmT))
  LinMod466209[i,1] <- x
  LinMod466209[i,2] <- try(anova(lmm2T, lmmT)[8][2,1])
  #print(try(anova(lmm2T, lmmT)[8][2,1]))
}


write.csv(LinMod466209,"ServerRunModel_77SamplesCD_inflamedvsNoninflamed.csv")

save(LinMod466209,file = "ServerRunModel_77SamplesCD_inflamedvsNoninflamed.RData")

# calculate adjusted p values
colnames(LinMod466209) <- c("CpG", "p_value")
LinMod466209$FDR_pval <- p.adjust(LinMod466209$p_value, "fdr")
LinMod466209$bonferroni_pval <- p.adjust(LinMod466209$p_value, "bonferroni")

write.csv(LinMod466209,"ServerRunModel_77Samples_All_fdr_bon_IvsNI.csv")

save(LinMod466209,file = "ServerRunModel_77Samples_All_fdr_bon_IvsNI.RData")

sink("ServerRunModel_77Samples_All_SessionInfo_IvsNI.txt")
sessionInfo()
dim(LinMod466209)
sink()



# get significant hits
LinTopHits_IvsNI <- LinMod466209[order(LinMod466209$FDR_pval),]

head(LinMod466209[order(LinMod466209$FDR_pval),],25)

LinTopHits_IvsNI_sig <- LinTopHits_IvsNI[LinTopHits_IvsNI$FDR_pval < 0.05,]

save(LinTopHits_IvsNI_sig,file = "LinTopHits_IvsNI_sig_FDR_BON.RData")



############################################ pc model


LinMod466209_IvsNI <- LinMod466209


load("LinTopHits_IvsNI_sig_FDR_BON.RData")

load("100_samples_updated_metadata.RData")


#load("pcaResandpc1to10plusMeta.Rdata")
# this is 100 samples, generate with 77 instead

metaAll_temp <- meta100
qqDF_temp <- qq.ready

topHits1 <- LinTopHits_IvsNI_sig

##  I vs NI

keepCGIvsNI <- topHits1$CpG
qqDF_temp1 <- qqDF_temp[keepCGIvsNI,]
dim(qqDF_temp1)

metaAll_temp1 <- metaAll_temp[metaAll_temp$Condition == "CD",]
qqDF_temp2 <- qqDF_temp1[,metaAll_temp1$SampleID.1]

dim(metaAll_temp1)
dim(qqDF_temp2)

pcareresults=prcomp(t(qqDF_temp2))

PCA_res <- pcareresults$x
PCA_res <- as.data.frame(PCA_res)
PCA_res$SampleID.1 <- rownames(PCA_res)


metaAll_temp2 <- merge(metaAll_temp, PCA_res, by="SampleID.1")


class(metaAll_temp2$Status)
table(metaAll_temp2$Status)


class(metaAll_temp2$Patient_ID)

metaAll_temp2Keep <- metaAll_temp2
library(lme4)


PC_ps10_IvsNI <- data.frame(matrix(, nrow=169138, ncol=2))

for (i in (1:169138)){
  #print(rownames(qqDF_temp1)[i])
  x <- rownames(qqDF_temp2)[i]
  metaAll_tempX <- metaAll_temp2
  metaAll_tempX$temp <- t(qqDF_temp2[i,])[1,]
  #print(t(qqDF_temp1[i,]))
  lmmT <- lmer(temp ~ Status + PC1 + PC2 + PC3 + PC4 + PC5 + PC6+ PC7 + PC8 + PC9 + PC10 + (1 | Patient_ID), data = metaAll_tempX,
               REML = FALSE)
  lmm2T <- lmer(temp ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6+ PC7 + PC8 + PC9 + PC10 + (1 | Patient_ID), data = metaAll_tempX,
                REML = FALSE)
  #print(summary(lmm))
  #print(anova(lmm2T, lmmT))
  PC_ps10_IvsNI[i,1] <- x
  PC_ps10_IvsNI[i,2] <- anova(lmm2T, lmmT)[8][2,1]
}


colnames(PC_ps10_IvsNI) <- c("CpG", "PC10_pValue")



save(PC_ps10_IvsNI, file ="PC10_linMod_IvsNI.Rdata")

# WARNINGS/ERRORS
#boundary (singular) fit: see ?isSingular
#Warning messages:
# 1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                  Model failed to converge with max|grad| = 0.00233573 (tol = 0.002, component 1)







library(dplyr)
library(readr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")

query_cpgs <- LinTopHits_IvsNI_sig$CpG


anno       <- IlluminaHumanMethylation450kanno.ilmn12.hg19 %>% 
  getAnnotation %>% 
  as.data.frame %>% 
  dplyr::slice(match(query_cpgs, Name))

dim(anno)
head(anno)
colnames(anno)

LinTopHits_IvsNI_sig$Name <- LinTopHits_IvsNI_sig$CpG
LinTopHits_IvsNI_anno <- merge(LinTopHits_IvsNI_sig, anno, by="Name")

save(LinTopHits_IvsNI_anno,file = "LinTopHits_IvsNI_anno_FDR_BON.RData")



