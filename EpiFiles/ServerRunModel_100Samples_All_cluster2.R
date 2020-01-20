# load in normalised betas -> qq.ready
# loading in metadata for 100 CD samples -> meta100


load("100_samples_updated_metadata.RData")
library(lme4)

LinMod466209 <- data.frame(matrix(, nrow=466209, ncol=2))

metaAll_temp <- meta100
qqDF_temp1 <- qq.ready

class(metaAll_temp$cluster)
table(metaAll_temp$cluster)
metaAll_temp$cluster[metaAll_temp$cluster != "2"] <- "Other"
metaAll_temp$cluster <- as.factor(metaAll_temp$cluster)
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)

# Run linear mixed model
for (i in (1:466209)){
  lmmT <- NULL
  lmm2T <- NULL
  print(rownames(qqDF_temp1)[i])
  x <- rownames(qqDF_temp1)[i]
  metaAll_temp1 <- metaAll_temp
  metaAll_temp1$temp <- t(qqDF_temp1[i,])[1,]
  #print(t(qqDF_temp1[i,]))
  lmmT <- try(lmer(temp ~ cluster + Sex + Condition_Status+  Smoker + Age + Sentrix_ID + Sentrix_Position  + Biopsy_location +Hospital+ (1 | Patient_ID), data = metaAll_temp1,
                   REML = FALSE))
  lmm2T <- try(lmer(temp ~ Sex + Condition_Status+ Smoker + Age + Sentrix_ID + Sentrix_Position  + Biopsy_location + Hospital+(1 | Patient_ID), data = metaAll_temp1,
                    REML = FALSE))
  #print(summary(lmm))
  #print(anova(lmm2T, lmmT))
  LinMod466209[i,1] <- x
  LinMod466209[i,2] <- try(anova(lmm2T, lmmT)[8][2,1])
  #print(try(anova(lmm2T, lmmT)[8][2,1]))
}


write.csv(LinMod466209,"ServerRunModel_100Samples_All_cluster2.csv")

save(LinMod466209,file = "ServerRunModel_100Samples_All_cluster2.RData")

# calculate adjusted p values
colnames(LinMod466209) <- c("CpG", "p_value")
LinMod466209$FDR_pval <- p.adjust(LinMod466209$p_value, "fdr")
LinMod466209$bonferroni_pval <- p.adjust(LinMod466209$p_value, "bonferroni")

write.csv(LinMod466209,"ServerRunModel_100Samples_All_fdr_pval_cluster2.csv")

save(LinMod466209,file = "ServerRunModel_100Samples_All_fdr_pval_cluster2.RData")

sink("ServerRunModel_100Samples_All_SessionInfo_cluster2.txt")
sessionInfo()
dim(LinMod466209)
sink()


library(dplyr)
library(readr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")


# get significant hits
LinTopHits_2 <- LinMod466209[order(LinMod466209$FDR_pval),]

head(LinMod466209[order(LinMod466209$FDR_pval),],814)

LinTopHits_2_sig <- LinTopHits_2[1:813,]
query_cpgs <- LinTopHits_2_sig$CpG

#get illumina annotation for sig hits
anno       <- IlluminaHumanMethylation450kanno.ilmn12.hg19 %>% 
  getAnnotation %>% 
  as.data.frame %>% 
  dplyr::slice(match(query_cpgs, Name))

dim(anno)
head(anno)
colnames(anno)

LinTopHits_2_sig$Name <- LinTopHits_2_sig$CpG
LinTopHits_2_anno <- merge(LinTopHits_2_sig, anno, by="Name")

save(LinTopHits_2_anno,file = "LinTopHits_2_anno_FDR_BON.RData")
write.csv(LinTopHits_2_anno,file = "LinTopHits_2_anno_FDR_BON.csv")

