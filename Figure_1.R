library(ggplot2)
library(phyloseq)
library(ape)
library(plyr)
library(gridExtra)

### Reading in Data

mapping = readRDS("SIRG_DADA2_mapping.RDS")
raw_cnts = readRDS("SIRG_DADA2_nc.nh.counts.RDS")
tree = read.tree("SIRG_biopsy_DADA2.final.nohuman.phylo.tree")
mothur = read.table("SIRG_biopsy_DADA2.mothur.txt",header=TRUE,sep="\t",row.names = 1)
spingo = read.table("SIRG_biopsy_DADA2.final.spingo.7.txt",header=TRUE,sep="\t",row.names=1)


#Filter counts for features present in at least 5% of samples
cnts_filt10 = raw_cnts[apply(raw_cnts>0,1,sum)>=round(ncol(raw_cnts)*0.05),]
map_data = sample_data(mapping) #For PhyloSeq


### Generate Relative Abundance data

cnts_filt10_prop = prop.table(as.matrix(cnts_filt10),2)
cnts_filt10_prop = cnts_filt10_prop * 100

### Getting sample with most reads per patient 
sizes = sort(apply(raw_cnts,2,sum),decreasing=TRUE)
toUnique = data.frame(sizes,"PatientID" = as.character(mapping[names(sizes),"PatientID"]))
toUnique = toUnique[!duplicated(toUnique$PatientID),]
bmap = mapping[rownames(toUnique),]
HTb = rownames(bmap[bmap$Condition=="Healthy",])


### Generating Alpha diversity

OTU_raw = otu_table(raw_cnts,taxa_are_rows = TRUE)
physeq_raw = phyloseq(OTU_raw,map_data,tree)
alpha_div = estimate_richness(physeq_raw)

### Generating Beta diversity

OTU = otu_table(cnts_filt10_prop,taxa_are_rows =TRUE)
physeq = phyloseq(OTU,map_data,tree)

### Bray Curtis distance in PhyloSeq
braycurtis_dist10 = phyloseq::distance(physeq,method="bray")
braycurtis_pc10_prop = pcoa(braycurtis_dist10)

### Making data frame for plotting in ggplot2 
BC_data.df = data.frame(mapping,
                        alpha_div,
                        "PC1" = braycurtis_pc10_prop$vectors[,1]*-1,
                        "PC2" = braycurtis_pc10_prop$vectors[,2]*-1)


HTs = rownames(BC_data.df[BC_data.df$Condition=="Healthy",])
HTn = HTs[!HTs %in% HTb]
BC_data.df2 = BC_data.df[!rownames(BC_data.df) %in% HTn,]
alpha_diversity = ggplot(BC_data.df2,aes(x=Description2,y=Shannon,fill=Condition))+geom_violin()+scale_fill_manual("Description2",values=c("black","#1e8bc3","#C30000","#ffb90f"))+theme_classic()+geom_dotplot(binaxis="y",stackdir = "center",dotsize=1,binwidth = 0.03,aes(fill="black"))

CDa = BC_data.df2[BC_data.df2$Description=="CD_Inflamed",]
CDi = BC_data.df2[BC_data.df2$Description=="CD_NonInflamed",]
UCa = BC_data.df2[BC_data.df2$Description=="UC_Inflamed",]
UCi = BC_data.df2[BC_data.df2$Description=="UC_Noninflamed",]
HTt = BC_data.df2[BC_data.df2$Description2=="Healthy",]



## PC1 correlations
PC1_BC_cors = apply(cnts_filt10_prop,1,function(x){
  cor(BC_data.df$PC1,x,method = "spearman")})

PC1_BC_pvalues = p.adjust(apply(cnts_filt10_prop,1,function(x){
  cor.test(BC_data.df$PC1,x,method = "spearman")$p.value}),method = "fdr")

PC1_pvalues = PC1_BC_pvalues[names(PC1_BC_cors)]
PC1_cors = data.frame("PC1_cors" = PC1_BC_cors,"PC1_pvalues" = PC1_pvalues)


## PC2 correlations
PC2_BC_cors = apply(cnts_filt10_prop,1,function(x){
  cor(BC_data.df$PC2,x,method = "spearman")})
PC2_BC_pvalues = p.adjust(apply(cnts_filt10_prop,1,function(x){
  cor.test(BC_data.df$PC2,x,method = "spearman")$p.value}),method = "fdr")

PC2_pvalues = PC2_BC_pvalues[names(PC2_BC_cors)]
PC2_cors = data.frame("PC2_cors" = PC2_BC_cors,"PC2_pvalues" = PC2_pvalues)



PC_cors = data.frame(mothur[rownames(PC1_cors),],spingo=spingo[rownames(PC1_cors),"Species"],
                     cbind(PC1_cors,PC2_cors))
PC_cors = PC_cors[(PC_cors$PC1_pvalues <0.05 | PC_cors$PC2_pvalues <0.05),]


### Correlation between PC1 and Shannon diversity 
cor(BC_data.df$PC1,alpha_div$Shannon,method = "spearman")
cor(BC_data.df$PC2,alpha_div$Shannon,method = "spearman")


### PCoA plot
pcoa_plot = ggplot(BC_data.df,aes(x = PC1, y = PC2,color=Condition))+
  geom_point(size=2,aes(color=Condition,shape=factor(mapping$Status)))+
  geom_point(size=2,aes(shape=factor(BC_data.df$Status)))+scale_shape_manual(values=c(1,17))+
  geom_line(data=BC_data.df,aes(x=PC1,y=PC2,group=mapping$PatientID))+
  stat_ellipse(aes(x = PC1,y=PC2,linetype=factor(mapping$Status),fill=factor(mapping$Condition)),geom="polygon",level=0.8,alpha=0.2)+
  scale_linetype_manual("Ellipse Linetype",values=c("Inflamed"=2,"NonInflamed"=1))+
  xlab(paste("PC1:",round(braycurtis_pc10_prop$values$Relative_eig[1]*100,2)," variation"))+
  ylab(paste("PC2:",round(braycurtis_pc10_prop$values$Relative_eig[2]*100,2)," variation"))+
  scale_fill_manual("Condition",values=c("#1e8bc3","#C30000","#ffb90f"))+
  scale_colour_manual("Condition",values=c("#1e8bc3","#C30000",	"#ffb90f"))+guides(color=FALSE,lty=FALSE,shape=guide_legend(title="Sample"),fill=guide_legend(title="Condition"))+theme(panel.background=element_blank())+
  ggtitle("Bray Curtis")

### Plotting correlation arrows
arrows = ggplot(BC_data.df,aes(x = PC1/max(abs(PC1)), y = PC2/max(abs(PC2)),color=Condition))+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = 0.06823481, yend = -0.7487282),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = 0.1386036, yend = 0.648119),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = -0.6881125, yend = -0.1166202),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = 0.5695726, yend = -0.01284271),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  geom_segment(color="green3",linetype=1,aes(x = 0, y = 0, xend = 0.6478033, yend = 0.06917217),size=1.5, arrow = arrow(length = unit(0.5, "cm")))+
  geom_text(x=0.06823481,y=-0.7887282,label="B. vulgatus",color="black",fontface="italic")+
  geom_text(x=-0.7081125,y=-0.1766202,label="Esch/Shig",color="black",fontface="italic")+
  geom_text(x=0.605726,y=-0.07284271,label="F. prausnitzii",color="black",fontface="italic")+
  geom_text(x=0.1086036,y=0.678119,label="B. dorei",color="black",fontface="italic")+
  geom_text(x=0.6278033,y=0.1417217,label="Shannon Diversity",color="green3",fontface="italic")+theme_classic()+
  geom_text(x=-0.4626701,y=0.1360894,label="Clostridium sp.",color="black",fontface="italic")+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = -0.4626701, yend = 0.1360894),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  geom_text(x=-0.1027414,y=0.2496268,label="P. asaccharolytica",color="black",fontface="italic")+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = -0.1027414, yend = 0.2496268),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  geom_text(x=-0.24486,y=0.17515,label="B. fragilis",color="black",fontface="italic")+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = -0.24486, yend = 0.17515),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
 
  geom_text(x=-0.4456838,y=0.004321966,label="R.gnavus",color="black",fontface="italic")+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = -0.4456838, yend = 0.004321966),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  
  geom_text(x=0.3454578,y=0.03604542,label="Anaerostipes hadrus",color="black",fontface="italic")+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = 0.3454578, yend = 0.03604542),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  geom_text(x=0.4544845,y=-0.09129027,label="Lachnospiraceae sp.",color="black")+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend =0.4544845 , yend = -0.09129027),size=.7, arrow = arrow(length = unit(0.5, "cm")))



### Subtractive PC1 values
#This was wrapped in a for loop for testing all components
PC1 = braycurtis_pc10_prop$vectors[,1]
# CD
inflamed = grep("sCD[0-9]+a",names(PC1),value=TRUE,perl=TRUE)
CD_PC1_diff = c() 
for (sample in inflamed){
  a = PC1[sample]
  i = PC1[gsub("a","i",sample)]
  bin_val = a > i
  CD_PC1_diff = c(CD_PC1_diff,a - i)
}

# UC
inflamed = grep("sUC[0-9]+a",names(PC1),value=TRUE,perl=TRUE)
UC_PC1_diff = c()
for (sample in inflamed){
  a = PC1[sample]
  i = PC1[gsub("a","i",sample)]
  bin_val = a > i
  UC_PC1_diff = c(UC_PC1_diff,a - i) 
}


#HT
inflamed = grep("sHT[0-9]+.1",names(PC1),value=TRUE,perl=TRUE)
HT_PC1_diff = c()
for (sample in inflamed){
  sample1 = gsub(".1","",sample,fixed=TRUE)
  sample2 = paste(sample1,as.character(sample(1:2,1)),sep=".")
  a = PC1[sample2]
  totest = substr(sample2,nchar(sample2), nchar(sample2))
  if (totest == 1){
    i = PC1[gsub("\\.1","\\.2",sample2)]}
  else {
    i = PC1[gsub("\\.2","\\.1",sample2)]
  }
  HT_PC1_diff = c(HT_PC1_diff,a - i)
}

CD_PC1 = data.frame(rep("CD",length(CD_PC1_diff)),CD_PC1_diff)
colnames(CD_PC1) = c("Condition","PC1")
UC_PC1 = data.frame(rep("UC",length(UC_PC1_diff)),UC_PC1_diff)
colnames(UC_PC1) = c("Condition","PC1")
HT_PC1 = data.frame(rep("HT",length(HT_PC1_diff)),HT_PC1_diff)
colnames(HT_PC1) = c("Condition","PC1")
BC_PC1_data.df = rbind(CD_PC1,HT_PC1,UC_PC1)

CDx = wilcox.test(CD_PC1$PC1)
UCx = wilcox.test(UC_PC1$PC1)

BC_PC1_diff= ggplot(BC_PC1_data.df,aes(x = Condition,y=PC1,fill=Condition))+
  geom_violin()+geom_boxplot(width=0.3,outlier.colour = NA)+scale_fill_manual("Condition",values=c("#1e8bc3","#C30000",	"#ffb90f"))+
  ggtitle("PC1 BC")+geom_hline(yintercept = 0)+theme_classic()





CDa_no = (CDa$Shannon)
CDa_no = CDa_no[!is.na(CDa_no)]

CDi_no = (CDi$Shannon)
CDi_no = CDi_no[!is.na(CDi_no)]

UCa_no = (UCa$Shannon)
UCa_no = UCa_no[!is.na(UCa_no)]

UCi_no = (UCi$Shannon)
UCi_no = UCi_no[!is.na(UCi_no)]

HTt_no = (HTt$Shannon)
HTt_no = HTt_no[!is.na(HTt_no)]

CDa_shan = data.frame(rep("CDa",length(CDa_no)),CDa_no)
colnames(CDa_shan) = c("Condition","shan")

CDi_shan = data.frame(rep("CDi",length(CDi_no)),CDi_no)
colnames(CDi_shan) = c("Condition","shan")

UCa_shan = data.frame(rep("UCa",length(UCa_no)),UCa_no)
colnames(UCa_shan) = c("Condition","shan")

UCi_shan = data.frame(rep("UCi",length(UCi_no)),UCi_no)
colnames(UCi_shan) = c("Condition","shan")

HTt_shan = data.frame(rep("HT",length(HTt_no)),HTt_no)
colnames(HTt_shan) = c("Condition","shan")


BC_shan_data.df = rbind(CDa_shan,CDi_shan,
                        UCa_shan,UCi_shan,HTt_shan)



alpha_violin = ggplot(BC_data.df2,aes(x=Description2,y=Shannon,fill=Condition))+geom_violin()+scale_fill_manual("Description2",values=c("black","#1e8bc3","#C30000","#ffb90f"))+theme_classic()+geom_dotplot(binaxis="y",stackdir = "center",dotsize=1,binwidth = 0.03,aes(fill="black"))

grid.arrange(pcoa_plot,arrows,alpha_violin,BC_PC1_diff)


