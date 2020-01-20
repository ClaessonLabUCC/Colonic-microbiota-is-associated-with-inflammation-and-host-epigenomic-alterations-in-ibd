load("../100_samples_updated_metadata.RData")
load("LinTopHits_3_anno_FDR_BON.Rdata") 
#OR
load("ServerRunModel_100Samples_All_fdr_pval_cluster3.RData")
LinTopHits_3_sig <- LinTopHits_3_anno

# calculate delta for cluster 3 (difference in means across CpG sites)
head(rownames(LinTopHits_3_sig))
rowKeeps3 <- LinTopHits_3_sig$CpG
qqDF_RowsK <- qqDF_temp1[rowKeeps3,]
dim(qqDF_RowsK)
head(rownames(qqDF_RowsK))
metaAll_temp1[metaAll_temp1$cluster == "3",]
colsmet3 <- metaAll_temp1[metaAll_temp1$cluster == "3",]
colsmet3K <- colsmet3$SampleID.1
qqDF_RowsK3 <- qqDF_RowsK[,colsmet3K]
dim(qqDF_RowsK3)
colsmetO <- metaAll_temp1[metaAll_temp1$cluster == "other",]
dim(colsmetO)
colsmetOK <- colsmetO$SampleID.1
qqDF_RowsKO <- qqDF_RowsK[,colsmetOK]
dim(qqDF_RowsKO)
qqDF_RowsK3_spare <- as.data.frame(qqDF_RowsK3)
qqDF_RowsK3_spare$delta <- rowMeans(qqDF_RowsK3) - rowMeans(qqDF_RowsKO)
qqDF_RowsK3_spare$CpG <- rownames(qqDF_RowsK3_spare)
head(qqDF_RowsK3_spare)
delta_cluster3_sig <- qqDF_RowsK3_spare[,16:17]



save(delta_cluster3_sig,file = "delta_cluster3_sig.RData")
write.csv(delta_cluster3_sig,file = "delta_cluster3_sig.csv")
