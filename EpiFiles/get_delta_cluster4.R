
load("../100_samples_updated_metadata.RData")
load("LinTopHits_4_anno_FDR_BON.Rdata") 
#OR
load("ServerRunModel_100Samples_All_fdr_pval_cluster4.RData")
LinTopHits_4_sig <- LinTopHits_4_anno



# calculate delta for cluster 4 (difference in means across CpG sites)
head(rownames(LinTopHits_4_sig))
rowKeeps4 <- LinTopHits_4_sig$CpG
qqDF_RowsK <- qqDF_temp1[rowKeeps4,]
dim(qqDF_RowsK)
head(rownames(qqDF_RowsK))
metaAll_temp1[metaAll_temp1$cluster == "4",]
colsmet4 <- metaAll_temp1[metaAll_temp1$cluster == "4",]
colsmet4K <- colsmet4$SampleID.1
qqDF_RowsK4 <- qqDF_RowsK[,colsmet4K]
dim(qqDF_RowsK4)
colsmetO <- metaAll_temp1[metaAll_temp1$cluster == "Other",]
dim(colsmetO)
colsmetOK <- colsmetO$SampleID.1
qqDF_RowsKO <- qqDF_RowsK[,colsmetOK]
dim(qqDF_RowsKO)
qqDF_RowsK4_spare <- as.data.frame(qqDF_RowsK4)
qqDF_RowsK4_spare$delta <- rowMeans(qqDF_RowsK4) - rowMeans(qqDF_RowsKO)
qqDF_RowsK4_spare$CpG <- rownames(qqDF_RowsK4_spare)
head(qqDF_RowsK4_spare)
dim(qqDF_RowsK4_spare)
delta_cluster4_sig <- qqDF_RowsK4_spare[,14:15]



save(delta_cluster4_sig,file = "delta_cluster4_sig.RData")
write.csv(delta_cluster4_sig,file = "delta_cluster4_sig.csv")
