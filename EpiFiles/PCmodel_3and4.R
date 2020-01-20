
# /home/x/R/bin/R
# PC model for cell composition
library(lme4)

#######################################
load("ServerRunModel_100Samples_All_fdr_pval_cluster3.RData")
ls()
LinMod466209_C3 <- LinMod466209
LinMod466209 <- NULL
ls()
load("ServerRunModel_100Samples_All_fdr_pval_cluster4.RData")
LinMod466209_C4 <- LinMod466209
ls()


LinMod466209_C3_O <- LinMod466209_C3[order(LinMod466209_C3$FDR_pval),]
LinMod466209_C4_O <- LinMod466209_C4[order(LinMod466209_C4$FDR_pval),]

####################################

load("LinTopHits_3_anno_FDR_BON.RData")
ls()
load("LinTopHits_4_anno_FDR_BON.RData")
ls()

LinTopHits_3_sig <- LinTopHits_3_anno[,2:5]
LinTopHits_4_sig <- LinTopHits_4_anno[,2:5]


topHits3 <- LinTopHits_3_sig
topHits4 <- LinTopHits_4_sig


load("../100_samples_updated_metadata.RData")
load("../pcaResandpc1to10plusMeta.Rdata")


metaAll_temp <- meta100
qqDF_temp <- qq.ready


completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

## cluster 3 vs all

keepCG3 <- topHits3$CpG
qqDF_temp3 <- qqDF_temp[keepCG3,]
dim(qqDF_temp3)

PCA_res <- pcareresults$x
PCA_res <- as.data.frame(PCA_res)
PCA_res$SampleID.1 <- rownames(PCA_res)


metaAll_temp3 <- merge(metaAll_temp, PCA_res, by="SampleID.1")

class(metaAll_temp3$cluster)
table(metaAll_temp3$cluster)
metaAll_temp3$cluster[metaAll_temp3$cluster != 3] <- "Other"
metaAll_temp3$cluster <- as.factor(metaAll_temp3$cluster)
class(metaAll_temp3$cluster)
table(metaAll_temp3$cluster)

class(metaAll_temp3$Patient_ID)

metaAll_temp3Keep <- metaAll_temp3
library(lme4)

PC_ps10_C3 <- data.frame(matrix(, nrow=280, ncol=2))

for (i in (1:280)){
  #print(rownames(qqDF_temp3)[i])
  x <- rownames(qqDF_temp3)[i]
  metaAll_tempX <- metaAll_temp3
  metaAll_tempX$temp <- t(qqDF_temp3[i,])[1,]
  #print(t(qqDF_temp3[i,]))
  lmmT <- lmer(temp ~ cluster + PC1 + PC2 + PC3 + PC4 + PC5 + PC6+ PC7 + PC8 + PC9 + PC10 + (1 | Patient_ID), data = metaAll_tempX,
               REML = FALSE)
  lmm2T <- lmer(temp ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6+ PC7 + PC8 + PC9 + PC10 + (1 | Patient_ID), data = metaAll_tempX,
                REML = FALSE)
  #print(summary(lmm))
  #print(anova(lmm2T, lmmT))
  PC_ps10_C3[i,1] <- x
  PC_ps10_C3[i,2] <- anova(lmm2T, lmmT)[8][2,1]
}


colnames(PC_ps10_C3) <- c("CpG", "PC10_pValue_C3")

# WARNINGS/ERRORS
#boundary (singular) fit: see ?isSingular
#Warning messages:
 # 1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
  #                  Model failed to converge with max|grad| = 0.00233573 (tol = 0.002, component 1)
                  

write.csv(PC_ps10_C3, "PC10_linMod_3Vall.csv", row.names = F)

C3_table <- merge(topHits3, PC_ps10_C3, by ="CpG")


write.csv(C3_table, "PCs_C3_table.csv", row.names = F)

#nullMod <- lmer(qqnorm(meth) ~ cluster + PC1 + PC2 + . + PC5/10 + random-effects) 
#vs null (no cluster in the model); anova(null, full



#####################################################################

## cluster 4 vs all
keepCG4 <- topHits4$CpG
qqDF_temp4 <- qqDF_temp[keepCG4,]
dim(qqDF_temp4)

PCA_res <- pcareresults$x
PCA_res <- as.data.frame(PCA_res)
PCA_res$SampleID.1 <- rownames(PCA_res)


metaAll_temp4 <- merge(metaAll_temp, PCA_res, by="SampleID.1")

class(metaAll_temp4$cluster)
table(metaAll_temp4$cluster)
metaAll_temp4$cluster[metaAll_temp4$cluster != 4] <- "Other"
metaAll_temp4$cluster <- as.factor(metaAll_temp4$cluster)
class(metaAll_temp4$cluster)
table(metaAll_temp4$cluster)

class(metaAll_temp4$Patient_ID)

metaAll_temp4Keep <- metaAll_temp4
library(lme4)


PC_ps10_C4 <- data.frame(matrix(, nrow=813, ncol=2))

for (i in (1:813)){
  x <- rownames(qqDF_temp4)[i]
  metaAll_tempX <- metaAll_temp4
  metaAll_tempX$temp <- t(qqDF_temp4[i,])[1,]
  #print(t(qqDF_temp4[i,]))
  lmmT <- lmer(temp ~ cluster + PC1 + PC2 + PC3 + PC4 + PC5 + PC6+ PC7 + PC8 + PC9 + PC10 + (1 | Patient_ID), data = metaAll_tempX,
               REML = FALSE)
  lmm2T <- lmer(temp ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6+ PC7 + PC8 + PC9 + PC10 + (1 | Patient_ID), data = metaAll_tempX,
                REML = FALSE)
  #print(summary(lmm))
  #print(anova(lmm2T, lmmT))
  PC_ps10_C4[i,1] <- x
  PC_ps10_C4[i,2] <- anova(lmm2T, lmmT)[8][2,1]
}


#colnames(PC_ps6) <- c("CpG", "PC6_pValue")
colnames(PC_ps10_C4) <- c("CpG", "PC10_pValue_C4")

# WARNINGS/ERRORS
#boundary (singular) fit: see ?isSingular
#Warning messages:
# 1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                  Model failed to converge with max|grad| = 0.00233573 (tol = 0.002, component 1)


write.csv(PC_ps10_C4, "PC10_linMod_4Vall.csv", row.names = F)

C4_table <- merge(topHits4, PC_ps10_C4, by ="CpG")

write.csv(C4_table, "PCs_C4_table.csv", row.names = F)
