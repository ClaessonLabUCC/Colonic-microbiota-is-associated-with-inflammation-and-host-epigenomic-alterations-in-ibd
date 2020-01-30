require(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
library(data.table)

#inflamed versus noninflamed CD
pheno_data = read.csv("pheno_data_full.csv", header=TRUE)
pheno_data <- subset(pheno_data, pheno_data$ID %in% list.files("ballgown_CD"))
pheno_data$Status <- as.factor(pheno_data$Status)
setattr(pheno_data$Status,"levels",c("Inflamed", "noninflamed"))
pheno_data$ID == list.files("ballgown_CD") 

bg_chrX = ballgown(dataDir = "ballgown_CD/", samplePattern = "s", pData=pheno_data)

####alternate to paper

bg_filt = subset(bg_chrX,"rowVars(texpr(bg_chrX)) >1",genomesubset=TRUE)

bg_table = texpr(bg_filt, 'all')
bg_gene_names = unique(bg_table[, 9:10])

gene_expression = as.data.frame(gexpr(bg_filt))


head(gene_expression)

data_colors = c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')

transcript_gene_table = indexes(bg_chrX)$t2g

length(row.names(transcript_gene_table))
length(unique(transcript_gene_table[,"g_id"]))

counts=table(transcript_gene_table[,"g_id"])
c_one = length(which(counts == 1))
c_more_than_one = length(which(counts > 1))
c_max = max(counts)
hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
legend("topright", legend_text, lty=NULL)

full_table <- texpr(bg_chrX , 'all')
hist(full_table$length, breaks=50, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")

min_nonzero=1
data_columns=c(1:6)

boxplot(log2(gene_expression[,1:30]+min_nonzero), col=data_colors, names=pheno_data$X, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for all 30 libraries")

x = gene_expression[,"FPKM.sCD08a"]
y = gene_expression[,"FPKM.sCD08i"]
plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (sCDo8a, Inflamed sample)", ylab="FPKM (SCD08i, NonInflamed sample)", main="Comparison of expression values for a pair of samples")
abline(a=0,b=1)
rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")

colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
smoothScatter(x=log2(x+min_nonzero), y=log2(y+min_nonzero), xlab="FPKM (sCDo8a, Inflamed sample)", ylab="FPKM (SCD08i, NonInflamed sample)", main="Comparison of expression values for a pair of samples", colramp=colors, nbin=200)


gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)

i = which(gene_expression[,"sum"] > 5)

r=cor(gene_expression[i,data_columns])
r 

d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
par(mfrow=c(1,1))
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes) for all libraries", xlim=c(-0.15,0.15), ylim=c(-0.15,0.15))
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], pheno_data$ID, col=data_colors)

results_genes = stattest(bg_filt, feature="gene", covariate="Status", meas="FPKM", getFC = TRUE)
results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))

sig=which(results_genes$pval<0.05)
results_genes[,"de"] = log2(results_genes[,"fc"])
hist(results_genes[sig,"de"], breaks=50, col="seagreen", xlab="log2(Fold change) Cluster 1&2 versus 3 ", main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topleft", "Fold-change > 4", lwd=2, lty=2)

group1 <- subset(pheno_data , pheno_data$Status== "noninflamed")
group2 <- subset(pheno_data , pheno_data$Status == "Inflamed")

g1 <- as.character(group1$ID)
g2 <- as.character(group2$ID)

gene_expression[,"Other Clusters"]=apply(gene_expression[,g1], 1, mean)
gene_expression[,"Cluster 5"]=apply(gene_expression[,g2], 1, mean)
x=log2(gene_expression[,"Other Clusters"]+min_nonzero)
y=log2(gene_expression[,"Cluster 5"]+min_nonzero)
plot(x=x, y=y, pch=16, cex=0.25, xlab="Other Clusters (log2)", ylab="Cluster 5 (log2)", main="Other Clusters & 3 FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

write.table(results_genes, file="inflam_noninflam_CD.txt", sep="\t", row.names=FALSE, quote=FALSE)
sigpi = which(results_genes[,"pval"]<0.05)
sigp = results_genes[sigpi,]
#sigde = which(abs(sigp[,"de"]) >= 2)
#sig_tn_de = sigp[sigde,]

o = order(sigp[,"qval"], -abs(sigp[,"de"]), decreasing=FALSE)
output2 = sigp[o,c("gene_name","id","fc","pval","qval","de")]
output2 <- subset(output2, output2$pval < 0.05)
output2 <- subset(output2, output2$gene_name != ".")
write.csv(output2, file="inflam_noninflam_CD.csv", row.names=FALSE, quote=FALSE)
#View selected columns of the first 25 lines of output
output2[1:25,c(1,4,5)]

#inflamed versus noninflamed UC

pheno_data = read.csv("pheno_data_full - Copy - Copy.csv", header=TRUE)
pheno_data <- subset(pheno_data, pheno_data$ID %in% list.files("ballgown_UC"))
pheno_data$Status <- as.factor(pheno_data$Status)
setattr(pheno_data$Status,"levels",c("Inflamed", "noninflamed"))
pheno_data$ID == list.files("ballgown_UC") 

bg_chrX = ballgown(dataDir = "ballgown_UC/", samplePattern = "s", pData=pheno_data)

####alternate to paper

bg_filt = subset(bg_chrX,"rowVars(texpr(bg_chrX)) >1",genomesubset=TRUE)

bg_table = texpr(bg_filt, 'all')
bg_gene_names = unique(bg_table[, 9:10])

gene_expression = as.data.frame(gexpr(bg_filt))


head(gene_expression)

data_colors = c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')

transcript_gene_table = indexes(bg_chrX)$t2g

length(row.names(transcript_gene_table))
length(unique(transcript_gene_table[,"g_id"]))

counts=table(transcript_gene_table[,"g_id"])
c_one = length(which(counts == 1))
c_more_than_one = length(which(counts > 1))
c_max = max(counts)
hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
legend("topright", legend_text, lty=NULL)

full_table <- texpr(bg_chrX , 'all')
hist(full_table$length, breaks=50, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")

min_nonzero=1
data_columns=c(1:6)

boxplot(log2(gene_expression[,1:30]+min_nonzero), col=data_colors, names=pheno_data$X, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for all 30 libraries")

x = gene_expression[,"FPKM.sCD08a"]
y = gene_expression[,"FPKM.sCD08i"]
plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (sCDo8a, Inflamed sample)", ylab="FPKM (SCD08i, NonInflamed sample)", main="Comparison of expression values for a pair of samples")
abline(a=0,b=1)
rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")

colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
smoothScatter(x=log2(x+min_nonzero), y=log2(y+min_nonzero), xlab="FPKM (sCDo8a, Inflamed sample)", ylab="FPKM (SCD08i, NonInflamed sample)", main="Comparison of expression values for a pair of samples", colramp=colors, nbin=200)


gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)

i = which(gene_expression[,"sum"] > 5)

r=cor(gene_expression[i,data_columns])
r 

d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
par(mfrow=c(1,1))
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes) for all libraries", xlim=c(-0.15,0.15), ylim=c(-0.15,0.15))
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], pheno_data$ID, col=data_colors)

results_genes = stattest(bg_filt, feature="gene", covariate="Status", meas="FPKM", getFC = TRUE)
results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))

sig=which(results_genes$pval<0.05)
results_genes[,"de"] = log2(results_genes[,"fc"])
hist(results_genes[sig,"de"], breaks=50, col="seagreen", xlab="log2(Fold change) Cluster 1&2 versus 3 ", main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topleft", "Fold-change > 4", lwd=2, lty=2)

group1 <- subset(pheno_data , pheno_data$Status== "noninflamed")
group2 <- subset(pheno_data , pheno_data$Status == "Inflamed")

g1 <- as.character(group1$ID)
g2 <- as.character(group2$ID)

gene_expression[,"Other Clusters"]=apply(gene_expression[,g1], 1, mean)
gene_expression[,"Cluster 5"]=apply(gene_expression[,g2], 1, mean)
x=log2(gene_expression[,"Other Clusters"]+min_nonzero)
y=log2(gene_expression[,"Cluster 5"]+min_nonzero)
plot(x=x, y=y, pch=16, cex=0.25, xlab="Other Clusters (log2)", ylab="Cluster 5 (log2)", main="Other Clusters & 3 FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

write.table(results_genes, file="inflam_noninflam_UC.txt", sep="\t", row.names=FALSE, quote=FALSE)
sigpi = which(results_genes[,"pval"]<0.05)
sigp = results_genes[sigpi,]
#sigde = which(abs(sigp[,"de"]) >= 2)
#sig_tn_de = sigp[sigde,]

o = order(sigp[,"qval"], -abs(sigp[,"de"]), decreasing=FALSE)
output2 = sigp[o,c("gene_name","id","fc","pval","qval","de")]
output2 <- subset(output2, output2$pval < 0.05)
output2 <- subset(output2, output2$gene_name != ".")
write.csv(output2, file="inflam_noninflam_UC.csv", row.names=FALSE, quote=FALSE)
#View selected columns of the first 25 lines of output
output2[1:25,c(1,4,5)]


#healthy versus Disease CD
pheno_data = read.csv("pheno_data_full - Copy - Copy.csv", header=TRUE)
pheno_data <- subset(pheno_data, pheno_data$ID %in% list.files("ballgown_disease_CD"))
pheno_data$Condition <- as.factor(pheno_data$Condition)
setattr(pheno_data$Status,"levels",c("CD", "Healthy"))
pheno_data$ID == list.files("ballgown_disease_CD") 

bg_chrX = ballgown(dataDir = "ballgown_disease_CD/", samplePattern = "s", pData=pheno_data)


bg_filt = subset(bg_chrX,"rowVars(texpr(bg_chrX)) >1",genomesubset=TRUE)

bg_table = texpr(bg_filt, 'all')
bg_gene_names = unique(bg_table[, 9:10])

gene_expression = as.data.frame(gexpr(bg_filt))


head(gene_expression)

data_colors = c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')

transcript_gene_table = indexes(bg_chrX)$t2g

length(row.names(transcript_gene_table))
length(unique(transcript_gene_table[,"g_id"]))

counts=table(transcript_gene_table[,"g_id"])
c_one = length(which(counts == 1))
c_more_than_one = length(which(counts > 1))
c_max = max(counts)
hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
legend("topright", legend_text, lty=NULL)

full_table <- texpr(bg_chrX , 'all')
hist(full_table$length, breaks=50, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")

min_nonzero=1
data_columns=c(1:6)

boxplot(log2(gene_expression[,1:30]+min_nonzero), col=data_colors, names=pheno_data$X, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for all 30 libraries")

x = gene_expression[,"FPKM.sCD08a"]
y = gene_expression[,"FPKM.sCD08i"]
plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (sCDo8a, Inflamed sample)", ylab="FPKM (SCD08i, NonInflamed sample)", main="Comparison of expression values for a pair of samples")
abline(a=0,b=1)
rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")

colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
smoothScatter(x=log2(x+min_nonzero), y=log2(y+min_nonzero), xlab="FPKM (sCDo8a, Inflamed sample)", ylab="FPKM (SCD08i, NonInflamed sample)", main="Comparison of expression values for a pair of samples", colramp=colors, nbin=200)

gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)

i = which(gene_expression[,"sum"] > 5)

r=cor(gene_expression[i,data_columns])
r 

d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
par(mfrow=c(1,1))
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes) for all libraries", xlim=c(-0.15,0.15), ylim=c(-0.15,0.15))
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], pheno_data$ID, col=data_colors)

results_genes = stattest(bg_filt, feature="gene", covariate="Condition", meas="FPKM", getFC = TRUE)
results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))

sig=which(results_genes$pval<0.05)
results_genes[,"de"] = log2(results_genes[,"fc"])
hist(results_genes[sig,"de"], breaks=50, col="seagreen", xlab="log2(Fold change) Cluster 1&2 versus 3 ", main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topleft", "Fold-change > 4", lwd=2, lty=2)

group1 <- subset(pheno_data , pheno_data$Condition== "CD")
group2 <- subset(pheno_data , pheno_data$Condition == "Healthy")

g1 <- as.character(group1$ID)
g2 <- as.character(group2$ID)

gene_expression[,"Other Clusters"]=apply(gene_expression[,g1], 1, mean)
gene_expression[,"Cluster 5"]=apply(gene_expression[,g2], 1, mean)
x=log2(gene_expression[,"Other Clusters"]+min_nonzero)
y=log2(gene_expression[,"Cluster 5"]+min_nonzero)
plot(x=x, y=y, pch=16, cex=0.25, xlab="Other Clusters (log2)", ylab="Cluster 5 (log2)", main="Other Clusters & 3 FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

write.table(results_genes, file="Healthyvsdisease_CD.txt", sep="\t", row.names=FALSE, quote=FALSE)
sigpi = which(results_genes[,"pval"]<0.05)
sigp = results_genes[sigpi,]
#sigde = which(abs(sigp[,"de"]) >= 2)
#sig_tn_de = sigp[sigde,]

o = order(sigp[,"qval"], -abs(sigp[,"de"]), decreasing=FALSE)
output2 = sigp[o,c("gene_name","id","fc","pval","qval","de")]
output2 <- subset(output2, output2$pval < 0.01)
output2 <- subset(output2, output2$gene_name != ".")
write.csv(output2, file="Healthyvsdisease_CD.csv", row.names=FALSE, quote=FALSE)
#View selected columns of the first 25 lines of output
output2[1:25,c(1,4,5)]

#####Healthy versus disease UC

pheno_data = read.csv("pheno_data_full - Copy - Copy.csv", header=TRUE)
pheno_data <- subset(pheno_data, pheno_data$ID %in% list.files("ballgown_disease_UC"))
pheno_data$Condition <- as.factor(pheno_data$Condition)
setattr(pheno_data$Status,"levels",c("UC", "Healthy"))
pheno_data$ID == list.files("ballgown_disease_UC") 

bg_chrX = ballgown(dataDir = "ballgown_disease_UC/", samplePattern = "s", pData=pheno_data)

####alternate to paper

bg_filt = subset(bg_chrX,"rowVars(texpr(bg_chrX)) >1",genomesubset=TRUE)

bg_table = texpr(bg_filt, 'all')
bg_gene_names = unique(bg_table[, 9:10])

gene_expression = as.data.frame(gexpr(bg_filt))


head(gene_expression)

data_colors = c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')

transcript_gene_table = indexes(bg_chrX)$t2g

length(row.names(transcript_gene_table))
length(unique(transcript_gene_table[,"g_id"]))

counts=table(transcript_gene_table[,"g_id"])
c_one = length(which(counts == 1))
c_more_than_one = length(which(counts > 1))
c_max = max(counts)
hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
legend("topright", legend_text, lty=NULL)

full_table <- texpr(bg_chrX , 'all')
hist(full_table$length, breaks=50, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")

min_nonzero=1
data_columns=c(1:6)

boxplot(log2(gene_expression[,1:30]+min_nonzero), col=data_colors, names=pheno_data$X, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for all 30 libraries")

x = gene_expression[,"FPKM.sCD08a"]
y = gene_expression[,"FPKM.sCD08i"]
plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (sCDo8a, Inflamed sample)", ylab="FPKM (SCD08i, NonInflamed sample)", main="Comparison of expression values for a pair of samples")
abline(a=0,b=1)
rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")

colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
smoothScatter(x=log2(x+min_nonzero), y=log2(y+min_nonzero), xlab="FPKM (sCDo8a, Inflamed sample)", ylab="FPKM (SCD08i, NonInflamed sample)", main="Comparison of expression values for a pair of samples", colramp=colors, nbin=200)

gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)

i = which(gene_expression[,"sum"] > 5)

r=cor(gene_expression[i,data_columns])
r 

d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
par(mfrow=c(1,1))
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes) for all libraries", xlim=c(-0.15,0.15), ylim=c(-0.15,0.15))
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], pheno_data$ID, col=data_colors)

results_genes = stattest(bg_filt, feature="gene", covariate="Condition", meas="FPKM", getFC = TRUE)
results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))

sig=which(results_genes$pval<0.05)
results_genes[,"de"] = log2(results_genes[,"fc"])
hist(results_genes[sig,"de"], breaks=50, col="seagreen", xlab="log2(Fold change) Cluster 1&2 versus 3 ", main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topleft", "Fold-change > 4", lwd=2, lty=2)

group1 <- subset(pheno_data , pheno_data$Condition== "UC")
group2 <- subset(pheno_data , pheno_data$Condition == "Healthy")

g1 <- as.character(group1$ID)
g2 <- as.character(group2$ID)

gene_expression[,"Other Clusters"]=apply(gene_expression[,g1], 1, mean)
gene_expression[,"Cluster 5"]=apply(gene_expression[,g2], 1, mean)
x=log2(gene_expression[,"Other Clusters"]+min_nonzero)
y=log2(gene_expression[,"Cluster 5"]+min_nonzero)
plot(x=x, y=y, pch=16, cex=0.25, xlab="Other Clusters (log2)", ylab="Cluster 5 (log2)", main="Other Clusters & 3 FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

write.table(results_genes, file="Healthyvsdisease_UC.txt", sep="\t", row.names=FALSE, quote=FALSE)
sigpi = which(results_genes[,"pval"]<0.05)
sigp = results_genes[sigpi,]
#sigde = which(abs(sigp[,"de"]) >= 2)
#sig_tn_de = sigp[sigde,]

o = order(sigp[,"qval"], -abs(sigp[,"de"]), decreasing=FALSE)
output2 = sigp[o,c("gene_name","id","fc","pval","qval","de")]
output2 <- subset(output2, output2$pval < 0.01)
output2 <- subset(output2, output2$gene_name != ".")
write.csv(output2, file="Healthyvsdisease_UC.csv", row.names=FALSE, quote=FALSE)
#View selected columns of the first 25 lines of output
output2[1:25,c(1,4,5)]

