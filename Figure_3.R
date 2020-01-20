library("ggplot2")
library("plyr")
library("ape")
library("ade4")
library("made4")
source("heatplot_generic.R")
library("reshape2")
library("gridExtra")
library("phyloseq")
library(dendextend)
library(phyloseq)
library(xlsx)

## GGplot Colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## 27 samples for mRNA
mRNA = c("sCD04a","sCD04i","sCD18a","sCD18i","sCD24a","sCD24i","sCD42a","sCD42i",
         "sCD55a","sCD55i","sHT50.1","sHT50.2","sHT53.1","sHT53.2","sUC03a","sUC03i",
         "sUC101a.2","sUC101i.1","sUC101i.2","sUC17a.T1","sUC17i.T1","sUC68a","sUC68i",
         "sUC75a","sUC75i","sUC93a","sUC93i")


## Reading in tables
mapping = readRDS("SIRG_DADA2_mapping.RDS")
raw_cnts = readRDS("SIRG_DADA2_nc.nh.counts.RDS")
cnts_filt10 = raw_cnts[apply(raw_cnts>0,1,sum)>=round(ncol(raw_cnts)*0.05),]
tree = read.tree("SIRG_biopsy_DADA2.final.nohuman.phylo.tree")
mothur = read.table("SIRG_biopsy_DADA2.mothur.txt",header=TRUE,sep="\t",row.names = 1)
spingo = read.table("SIRG_biopsy_DADA2.final.spingo.7.txt",header=TRUE,sep="\t",row.names=1)

## Relative Abundance tables

cnts_filt10_prop = prop.table(as.matrix(cnts_filt10),2)
cnts_filt10_prop = cnts_filt10_prop * 100

### RSVs present in at least 10% of samples
#### Unifrac distances in PhyloSeq
OTU = otu_table(cnts_filt10_prop,taxa_are_rows =TRUE)
map_data = sample_data(mapping) #For PhyloSeq
physeq = phyloseq(OTU,map_data,tree)
braycurtis_dist10 = phyloseq::distance(physeq,method="bray")
braycurtis_pc10_prop = pcoa(braycurtis_dist10)


#### Preparing Hierarchical clustering ####
fac = as.factor(mapping$Description2)

########################
##### Bray Curtis ######
########################


par(oma=c(3,3,3,5))
#pdf(file="Figure_3_raw_v2.pdf",paper="A4r",width=11,height=8.5,onefile=TRUE,useDingbats=FALSE)
BC_OTUtree = heatplot_editd(cnts_filt10_prop,dist_in = braycurtis_dist10,scale='none',
                           method='ward.D2',
                           classvec=fac,
                           classvecCol = c("#1e8bc3","#1e8bc3","#C30000","#ffb90f","#ffb90f"),
                           returnSampleTree = TRUE)

library(dendextend)
BC_OTUtree_order = order.dendrogram(BC_OTUtree)
BC_OTUtree_samples = rownames(mapping)[BC_OTUtree_order]
BC_tree_clusters = cutree(BC_OTUtree,k=10) #Cut for 10 clusters, number determined with DynanmicTreeCut R package
BC_clusters = data.frame(clusters=BC_tree_clusters)
BC_clusters$clusters = factor(BC_clusters$clusters,
                              levels=c("1","2","3","4","5","6","7","8","9","10"))

cluster_cols = c("#F8766D","#00BFC4","#A3A500","#E76BF3","#8B4513","#D89000","#00B0F6","#FF62BC","#00BF7D","#9590FF")

heatplot_editd(cnts_filt10_prop,dist_in = braycurtis_dist10,scale='none',
               method='ward.D2',
               classvec=BC_clusters$clusters,
               classvecCol = cluster_cols,
               main="BrayCurtis Heatplot at OTU ward D2 with clusters")




######################### Family Counts ##################################

treeorder = labels(BC_OTUtree)
family = data.frame(mothur[rownames(raw_cnts),"Family"],raw_cnts)
rownames(family) = NULL
colnames(family)[1] = "Family"
family_con = ddply(family,"Family",numcolwise(sum))
rownames(family_con) = family_con[,1]
family_con = family_con[,-1]
family_con_prop = prop.table(as.matrix(family_con),2)
family_con_prop = family_con_prop * 100
tobesummed = rownames(family_con_prop)[apply(family_con_prop,1,mean)<1]
family_prop_tree = data.frame(as.character(rownames(family_con_prop)),family_con_prop)
colnames(family_prop_tree)[1] = "Family"
family_prop_tree$Family = as.character(family_prop_tree$Family)
tobesummed = tobesummed[-9]
family_prop_tree[tobesummed,"Family"] = "Other"
family_con_prop_tree_summed = ddply(family_prop_tree,"Family",numcolwise(sum))
rownames(family_con_prop_tree_summed) = family_con_prop_tree_summed$Family
family_con_prop_tree_summed = family_con_prop_tree_summed[,treeorder]

######################### Phylum Counts ##################################

phylum = data.frame(mothur[rownames(raw_cnts),"Phylum"],raw_cnts)
rownames(phylum) = NULL
colnames(phylum)[1] = "Phylum"
phylum_con = ddply(phylum,"Phylum",numcolwise(sum))
rownames(phylum_con) = phylum_con[,1]
phylum_con = phylum_con[,-1]
phylum_con_prop = prop.table(as.matrix(phylum_con),2)
phylum_con_prop = phylum_con_prop * 100
sort(apply(phylum_con_prop,1,sum),decreasing=TRUE) #Which phylym is highest


################### Reorder family table ###############################

taxa_lu = unique(mothur[,c("Family","Phylum")])
taxa_lu = taxa_lu[!taxa_lu$Family=="unclassified",]
rownames(taxa_lu) = taxa_lu$Family

phyla_order = as.character(taxa_lu[rownames(family_con_prop_tree_summed),"Phylum"])
phyla_order_sub = gsub("Firmicutes",1,phyla_order)
phyla_order_sub = gsub("Bacteroidetes",2,phyla_order_sub)
phyla_order_sub = gsub("Proteobacteria",3,phyla_order_sub)
phyla_order_sub = gsub("Actinobacteria",4,phyla_order_sub)

phyla_order_sub[7] = "5" #Replace NA from indexing with Other with position
phyla_order_sub[12] = "6" #Replace NA from indexing with Other with position
as.numeric(phyla_order_sub)
family_phyla = data.frame(phyla_order_sub,family_con_prop_tree_summed)
family_phyla = family_phyla[order(family_phyla$phyla_order),]
family_phyla2 = family_phyla[,-1]
apply(family_phyla2,1,sum) #Manually specify index based on what families you want in what order
family_phyla3 = family_phyla2[c(3,4,2,1,5,7,6,8,9,10,11,12),]

family_phyla3["Family"] = rownames(family_phyla3)
data.df = melt(family_phyla3,id.vars = "Family")
data.df$Family = factor(data.df$Family,levels=rev(rownames(family_phyla3)))


ggplot(data.df, aes(x = variable, y = value, fill = Family)) + geom_bar(stat = "identity",width=1) +
  scale_fill_manual(values = rev(c("red2","darkred","firebrick1","tomato2",
                               "blue1","darkslategray2","dodgerblue3",
                               "forestgreen","darkolivegreen3","goldenrod",
                               "grey","black")))+theme_classic()


