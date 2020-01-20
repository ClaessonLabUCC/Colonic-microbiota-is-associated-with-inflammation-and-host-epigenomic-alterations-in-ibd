################################
## Adjust fdr across all clusters as requested by reviewer
# can be found in the supplmentary table
c1 <- read.csv("ServerRunModel_100Samples_All_cluster1.csv")
c2 <- read.csv("ServerRunModel_100Samples_All_cluster2.csv")
c3 <- read.csv("ServerRunModel_100Samples_All_cluster3.csv")
c4 <- read.csv("ServerRunModel_100Samples_All_cluster4.csv")

head(c1)
head(c2)
head(c3)
head(c4)

colnames(c1) <- c("Row","CpG", "p_value")
colnames(c2) <- c("Row","CpG", "p_value")
colnames(c3) <- c("Row","CpG", "p_value")
colnames(c4) <- c("Row","CpG", "p_value")

c1 <- c1[,c(2:3)]
c2 <- c2[,c(2:3)]
c3 <- c3[,c(2:3)]
c4 <- c4[,c(2:3)]

head(c1)
head(c2)
head(c3)
head(c4)

c_all <- rbind(c1,c2,c3,c4)
head(c_all)

c_all$FDR_pval <- p.adjust(c_all$p_value, "fdr")
dim(c_all)
c_allsub1 <- c_all[c_all$p_value <0.05,]
dim(c_allsub1)

# all top cpgs
allTop <- read.csv("topHitsFDRsigPCsig_C1234.csv")

dim(allTop)
head(allTop)

Keep <-(allTop$CpG)

c_all_subset <- c_allsub1[c_allsub1$CpG %in% Keep, ]
dim(c_all_subset)
c_all_subset <- c_allsub1[c_allsub1$CpG %in% Keep, ]
dim(c_all_subset)


#########################################################

c_all_subset_C1 <- c_all_subset
colnames(c_all_subset_C1) <- c("CpG", "p_value","FDR_ALL_C1")
c_all_subset_C1$p_value_cluster1 <- c_all_subset_C1$p_value
c_all_subset_C1$p_value1 <- c_all_subset_C1$p_value
c_all_subset_C1$p_value <- NULL
all_clust1 <- merge(allTop, c_all_subset_C1, by=c("CpG", "p_value_cluster1"), all.x = T)
dim(allTop)
dim(all_clust1)

View(all_clust1[,c(1,2,5,56:57)])

###############################

c_all_subset_C2 <- c_all_subset
colnames(c_all_subset_C2) <- c("CpG", "p_value","FDR_ALL_C2")
c_all_subset_C2$p_value_cluster2 <- c_all_subset_C2$p_value
c_all_subset_C2$p_value2 <- c_all_subset_C2$p_value
c_all_subset_C2$p_value <- NULL
all_clust2 <- merge(all_clust1, c_all_subset_C2, by=c("CpG", "p_value_cluster2"), all.x = T)
dim(all_clust1)
dim(all_clust2)

View(all_clust2[,c(1,2,8,56:59)])


#################################
c_all_subset_C3 <- c_all_subset
colnames(c_all_subset_C3) <- c("CpG", "p_value","FDR_ALL_C3")
c_all_subset_C3$p_value_cluster3 <- c_all_subset_C3$p_value
c_all_subset_C3$p_value3 <- c_all_subset_C3$p_value
c_all_subset_C3$p_value <- NULL
all_clust3 <- merge(all_clust2, c_all_subset_C3, by=c("CpG", "p_value_cluster3"), all.x = T)
dim(all_clust2)
dim(all_clust3)

View(all_clust3[,c(1,2,11,55:61)])

##################################################################
c_all_subset_C4 <- c_all_subset
colnames(c_all_subset_C4) <- c("CpG", "p_value","FDR_ALL_C4")
c_all_subset_C4$p_value_cluster4 <- c_all_subset_C4$p_value
c_all_subset_C4$p_value4 <- c_all_subset_C4$p_value
c_all_subset_C4$p_value <- NULL
all_clust4 <- merge(all_clust3, c_all_subset_C4, by=c("CpG", "p_value_cluster4"), all.x = T)
dim(allTop)
dim(all_clust4)

View(all_clust4[,c(1,2,14,55:63)])

colnames(all_clust4)

all_clst <- all_clust4[,c(1,2,3,4,5,8,10,12,14,56:63)]
View(all_clst)


