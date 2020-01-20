# PC model for cell composition
load("ServerRunModel_100Samples_All_fdr_pval_cluster1.RData")
LinMod466209_C1 <- LinMod466209
LinMod466209 <- NULL
load("ServerRunModel_100Samples_All_fdr_pval_cluster2.RData")
LinMod466209_C2 <- LinMod466209

load("topHits1and2.RData")

load("../100_samples_updated_metadata.RData")
load("../pcaResandpc1to10plusMeta.Rdata")


metaAll_temp <- meta100
qqDF_temp <- qq.ready


completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

## cluster 1 vs all

keepCG1 <- topHits1$CpG
qqDF_temp1 <- qqDF_temp[keepCG1,]
dim(qqDF_temp1)

PCA_res <- pcareresults$x
PCA_res <- as.data.frame(PCA_res)
PCA_res$SampleID.1 <- rownames(PCA_res)


metaAll_temp1 <- merge(metaAll_temp, PCA_res, by="SampleID.1")

class(metaAll_temp1$cluster)
table(metaAll_temp1$cluster)
metaAll_temp1$cluster[metaAll_temp1$cluster != 1] <- "Other"
metaAll_temp1$cluster <- as.factor(metaAll_temp1$cluster)
class(metaAll_temp1$cluster)
table(metaAll_temp1$cluster)

class(metaAll_temp1$Patient_ID)

metaAll_temp1Keep <- metaAll_temp1
library(lme4)


PC_ps10_C1 <- data.frame(matrix(, nrow=817, ncol=2))

# linear model using principal components 
for (i in (1:817)){
  #print(rownames(qqDF_temp1)[i])
  x <- rownames(qqDF_temp1)[i]
  metaAll_tempX <- metaAll_temp1
  metaAll_tempX$temp <- t(qqDF_temp1[i,])[1,]
  #print(t(qqDF_temp1[i,]))
  lmmT <- lmer(temp ~ cluster + PC1 + PC2 + PC3 + PC4 + PC5 + PC6+ PC7 + PC8 + PC9 + PC10 + (1 | Patient_ID), data = metaAll_tempX,
               REML = FALSE)
  lmm2T <- lmer(temp ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6+ PC7 + PC8 + PC9 + PC10 + (1 | Patient_ID), data = metaAll_tempX,
                REML = FALSE)
  #print(summary(lmm))
  #print(anova(lmm2T, lmmT))
  PC_ps10_C1[i,1] <- x
  PC_ps10_C1[i,2] <- anova(lmm2T, lmmT)[8][2,1]
}


#colnames(PC_ps6) <- c("CpG", "PC6_pValue")
colnames(PC_ps10_C1) <- c("CpG", "PC10_pValue")

# WARNINGS/ERRORS
#boundary (singular) fit: see ?isSingular
#Warning messages:
 # 1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
  #                  Model failed to converge with max|grad| = 0.00233573 (tol = 0.002, component 1)
                  

write.csv(PC_ps10_C1, "PC10_linMod_1Vall.csv", row.names = F)

C1_table <- merge(topHits1, PC_ps10_C1, by ="CpG")
#C1_table <- merge(C1_table, PC_ps6, by ="CpG")

write.csv(C1_table, "C1_table.csv", row.names = F)


#nullMod <- lmer(qqnorm(meth) ~ cluster + PC1 + PC2 + . + PC5/10 + random-effects) 
#vs null (no cluster in the model); anova(null, full



#####################################################################

## cluster 2 vs all
keepCG2 <- topHits2$CpG
qqDF_temp2 <- qqDF_temp[keepCG2,]
dim(qqDF_temp2)

PCA_res <- pcareresults$x
PCA_res <- as.data.frame(PCA_res)
PCA_res$SampleID.1 <- rownames(PCA_res)


metaAll_temp2 <- merge(metaAll_temp, PCA_res, by="SampleID.1")

class(metaAll_temp2$cluster)
table(metaAll_temp2$cluster)
metaAll_temp2$cluster[metaAll_temp2$cluster != 2] <- "Other"
metaAll_temp2$cluster <- as.factor(metaAll_temp2$cluster)
class(metaAll_temp2$cluster)
table(metaAll_temp2$cluster)

class(metaAll_temp2$Patient_ID)

metaAll_temp2Keep <- metaAll_temp2
library(lme4)


PC_ps10_C2 <- data.frame(matrix(, nrow=11423, ncol=2))

for (i in (1:11423)){
  x <- rownames(qqDF_temp2)[i]
  metaAll_tempX <- metaAll_temp2
  metaAll_tempX$temp <- t(qqDF_temp2[i,])[1,]
  #print(t(qqDF_temp2[i,]))
  lmmT <- lmer(temp ~ cluster + PC1 + PC2 + PC3 + PC4 + PC5 + PC6+ PC7 + PC8 + PC9 + PC10 + (1 | Patient_ID), data = metaAll_tempX,
               REML = FALSE)
  lmm2T <- lmer(temp ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6+ PC7 + PC8 + PC9 + PC10 + (1 | Patient_ID), data = metaAll_tempX,
                REML = FALSE)
  #print(summary(lmm))
  #print(anova(lmm2T, lmmT))
  PC_ps10_C2[i,1] <- x
  PC_ps10_C2[i,2] <- anova(lmm2T, lmmT)[8][2,1]
}


#colnames(PC_ps6) <- c("CpG", "PC6_pValue")
colnames(PC_ps10_C2) <- c("CpG", "PC10_pValue_C2")

# WARNINGS/ERRORS
#boundary (singular) fit: see ?isSingular
#Warning messages:
# 1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                  Model failed to converge with max|grad| = 0.00233573 (tol = 0.002, component 1)


write.csv(PC_ps10_C2, "PC10_linMod_2Vall.csv", row.names = F)

C2_table <- merge(topHits2, PC_ps10_C2, by ="CpG")

write.csv(C2_table, "PCs_C2_table.csv", row.names = F)
