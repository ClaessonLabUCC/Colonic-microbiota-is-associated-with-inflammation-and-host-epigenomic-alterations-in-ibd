
library(dplyr)
library(readr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")

load("ServerRunModel_100Samples_All_fdr_pval_cluster1.RData")
LinMod466209_C1 <- LinMod466209
LinMod466209 <- NULL
load("ServerRunModel_100Samples_All_fdr_pval_cluster2.RData")
LinMod466209_C2 <- LinMod466209


LinMod466209_C1_O <- LinMod466209_C1[order(LinMod466209_C1$FDR_pval),]
LinMod466209_C2_O <- LinMod466209_C2[order(LinMod466209_C2$FDR_pval),]


head(LinMod466209_C1[order(LinMod466209_C1$FDR_pval),],818)
LinTopHits_1_sig <- LinMod466209_C1_O[1:817,]

head(LinMod466209_C2_O,11425)
LinTopHits_2_sig <- LinMod466209_C2_O[1:11423,]

topHits1 <- LinTopHits_1_sig
topHits2 <- LinTopHits_2_sig

colnames(topHits1)[2] <- "p_value_cluster1"
colnames(topHits2)[2] <- "p_value_cluster2"
colnames(topHits1)[3] <- "FDR_pval_cluster1"
colnames(topHits2)[3] <- "FDR_pval_cluster2"
colnames(topHits1)[4] <- "bonferroni_pval_cluster1"
colnames(topHits2)[4] <- "bonferroni_pval_cluster2"

topHits1and2 <- merge(topHits1, topHits2, by ="CpG", all=T)
dim(topHits1and2)

query_cpgs <- topHits1and2$CpG

anno       <- IlluminaHumanMethylation450kanno.ilmn12.hg19 %>% 
  getAnnotation %>% 
  as.data.frame %>% 
  dplyr::slice(match(query_cpgs, Name))

dim(anno)
head(anno)
colnames(anno)

topHits1and2$Name <- topHits1and2$CpG
topHits1and2_anno <- merge(topHits1and2, anno, by="Name")

save(topHits1and2_anno,file = "topHits1and2_anno.RData")
write.csv(topHits1and2_anno,file = "topHits1and2_anno.csv")


save(topHits1, topHits2,file = "topHits1and2.RData")
