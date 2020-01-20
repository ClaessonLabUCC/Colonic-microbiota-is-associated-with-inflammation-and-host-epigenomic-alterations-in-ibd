
load("pcaResandpc1to10plusMeta.Rdata")
results=round(((pcareresults$sdev)^2 / sum(pcareresults$sdev^2))*100,digits=2)
labels=as.character(results[1:10])
labels=paste(labels,"%",sep="")
barplot(results[1:10], names=labels,ylim=c(0,25 ),col="yellow")


#ggplot2 for PC1 vs PC2, colouring for condition, shapes for status and lines to connect samples from same patient.
labx1=paste("PC1",labels[1],sep=" ")
laby1=paste("PC2",labels[2],sep=" ")
laby2=paste("PC4",labels[4],sep=" ")
lab3=paste("PC3",labels[3],sep=" ")
lab5=paste("PC5",labels[5],sep=" ")
lab6=paste("PC6",labels[6],sep=" ")
lab7=paste("PC7",labels[7],sep=" ")
lab8=paste("PC8",labels[8],sep=" ")
lab9=paste("PC9",labels[9],sep=" ")
lab10=paste("PC10",labels[10],sep=" ")

data.df$cs <- paste(data.df$Condition, data.df$Status, sep=" ")
data.df$Inflamation <- factor(data.df$cs, levels = c("HT Non-inflamed", "CD Non-inflamed", "CD Inflamed"))
data.df$PC1byMinus1 <- -1*data.df$PC1



ggplot(data.df,aes(x = PC1byMinus1, y = PC2,color=Condition))+
  geom_point(size=2,aes(color=Condition,shape=Status))+ geom_line(aes(group=Patient_ID))+
  labs(x=labx1,y=laby1 )+
  theme_bw() +
  stat_ellipse(aes(lty=Inflamation, fill=Condition), geom = "polygon", level = 0.8, alpha=0.2, size=0.4) +
  scale_linetype_manual(values = c(1,1,2)) +
  scale_fill_manual("Condition",values=c("#1e8bc3","#C30000","#ffb90f"))+
  scale_color_manual("Condition",values=c("#1e8bc3","#C30000","#ffb90f"))+guides(color=FALSE)+
  theme(legend.position = "none",panel.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


ggplot(data.df,aes(x = PC1byMinus1, y = PC6,color=Condition))+
  geom_point(size=2,aes(color=Condition,shape=Status))+ geom_line(aes(group=Patient_ID))+
  labs(x=labx1,y=lab6 )+
  theme_bw() +
  stat_ellipse(aes(lty=Inflamation, fill=Condition), geom = "polygon", level = 0.8, alpha=0.2, size=0.4) +
  scale_linetype_manual(values = c(1,1,2)) +
  scale_fill_manual("Condition",values=c("#1e8bc3","#C30000","#ffb90f"))+
  scale_color_manual("Condition",values=c("#1e8bc3","#C30000","#ffb90f"))+guides(color=FALSE)+
  theme(legend.position = "none",panel.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
