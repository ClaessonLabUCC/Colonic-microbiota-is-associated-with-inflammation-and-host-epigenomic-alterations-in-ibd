
load("../100_samples_updated_metadata.RData")
load("topHits1and2") 
#OR
load("ServerRunModel_100Samples_All_fdr_pval_cluster2.RData")
LinTopHits_2_sig <- topHits2

# calculate delta for cluster 2 (difference in means across CpG sites)
metaAll_temp <- meta100
qqDF_temp1 <- qq.ready


class(metaAll_temp$cluster)
table(metaAll_temp$cluster)
metaAll_temp$cluster[metaAll_temp$cluster != "2"] <- "Other"
metaAll_temp$cluster <- as.factor(metaAll_temp$cluster)
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)

metaAll_temp1 <- metaAll_temp

head(rownames(LinTopHits_2_sig))
rowKeeps2 <- LinTopHits_2_sig$CpG
qqDF_RowsK <- qqDF_temp1[rowKeeps2,]
dim(qqDF_RowsK)
head(rownames(qqDF_RowsK))
metaAll_temp1[metaAll_temp1$cluster == "2",]
colsmet2 <- metaAll_temp1[metaAll_temp1$cluster == "2",]
colsmet2K <- colsmet2$SampleID.1
qqDF_RowsK2 <- qqDF_RowsK[,colsmet2K]
dim(qqDF_RowsK2)
colsmetO <- metaAll_temp1[metaAll_temp1$cluster == "Other",]
dim(colsmetO)
colsmetOK <- colsmetO$SampleID.1
qqDF_RowsKO <- qqDF_RowsK[,colsmetOK]
dim(qqDF_RowsKO)
qqDF_RowsK2_spare <- as.data.frame(qqDF_RowsK2)
qqDF_RowsK2_spare$delta <- rowMeans(qqDF_RowsK2) - rowMeans(qqDF_RowsKO)
qqDF_RowsK2_spare$CpG <- rownames(qqDF_RowsK2_spare)
head(qqDF_RowsK2_spare)
delta_cluster2_sig <- qqDF_RowsK2_spare[,13:14]



save(delta_cluster2_sig,file = "delta_cluster2_sig.RData")
write.csv(delta_cluster2_sig,file = "delta_cluster2_sig.csv")
