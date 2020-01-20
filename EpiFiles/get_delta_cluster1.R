
load("../100_samples_updated_metadata.RData")
load("topHits1and2") 
#OR
load("ServerRunModel_100Samples_All_fdr_pval_cluster1.RData")
LinTopHits_1_sig <- topHits1

# calculate delta for cluster 1 (difference in means across CpG sites)

metaAll_temp <- meta100
qqDF_temp1 <- qq.ready


class(metaAll_temp$cluster)
table(metaAll_temp$cluster)
metaAll_temp$cluster[metaAll_temp$cluster != "1"] <- "Other"
metaAll_temp$cluster <- as.factor(metaAll_temp$cluster)
class(metaAll_temp$cluster)
table(metaAll_temp$cluster)

metaAll_temp1 <- metaAll_temp

head(rownames(LinTopHits_1_sig))
rowKeeps1 <- LinTopHits_1_sig$CpG
qqDF_RowsK <- qqDF_temp1[rowKeeps1,]
dim(qqDF_RowsK)
head(rownames(qqDF_RowsK))
metaAll_temp1[metaAll_temp1$cluster == "1",]
colsmet1 <- metaAll_temp1[metaAll_temp1$cluster == "1",]
colsmet1K <- colsmet1$SampleID.1
qqDF_RowsK1 <- qqDF_RowsK[,colsmet1K]
dim(qqDF_RowsK1)
colsmetO <- metaAll_temp1[metaAll_temp1$cluster == "Other",]
dim(colsmetO)
colsmetOK <- colsmetO$SampleID.1
qqDF_RowsKO <- qqDF_RowsK[,colsmetOK]
dim(qqDF_RowsKO)
qqDF_RowsK1_spare <- as.data.frame(qqDF_RowsK1)
qqDF_RowsK1_spare$delta <- rowMeans(qqDF_RowsK1) - rowMeans(qqDF_RowsKO)
qqDF_RowsK1_spare$CpG <- rownames(qqDF_RowsK1_spare)
head(qqDF_RowsK1_spare)
delta_cluster1_sig <- qqDF_RowsK1_spare[,12:13]



save(delta_cluster1_sig,file = "delta_cluster1_sig.RData")
write.csv(delta_cluster1_sig,file = "delta_cluster1_sig.csv")
