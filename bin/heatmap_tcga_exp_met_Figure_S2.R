rm(list = ls())
# library(limma)
# 
# # rt3=read.table("tcgaSymbol_LGG.txt",sep="\t",header=T,check.names=F)
# rt3=read.table("tcgaSymbol_GBM.txt",sep="\t",header=T,check.names=F)
# rt3=as.matrix(rt3)
# rownames(rt3)=rt3[,1]
# exp3=rt3[,2:ncol(rt3)]
# dimnames3=list(rownames(exp3),colnames(exp3))
# data3=matrix(as.numeric(as.matrix(exp3)),nrow=nrow(exp3),dimnames=dimnames3)
# data3=avereps(data3)
# data3=data3[rowMeans(data3)>0.5,]
# dim(data3)
# group=sapply(strsplit(colnames(data3),"\\-"), "[", 4)
# group=sapply(strsplit(group,""), "[", 1)
# group=gsub("2", "1", group)
# table(group)
# data3=data3[,group==0]
# colnames(data3)
# 
# # met <- read.table("NID2_TCGA-LGG.met.txt",sep="\t",header=T,check.names=F)
# met <- read.table("NID2_TCGA-GBM_450.met.txt", sep="\t",header=T,check.names=F)
# met=as.matrix(met)
# rownames(met)=met[,1]
# met=met[,2:ncol(met)]
# dimnames3=list(rownames(met),colnames(met))
# met=matrix(as.numeric(as.matrix(met)),nrow=nrow(met),dimnames=dimnames3)
# met=avereps(met)
# group=sapply(strsplit(colnames(met),"\\-"), "[", 4)
# group=sapply(strsplit(group,""), "[", 1)
# group=gsub("2", "1", group)
# table(group)
# met=met[,group==0]
# colnames(met)
# 
# met2 <- read.table("NID2_TCGA-GBM_27.met.txt", sep="\t",header=T,check.names=F)
# met2=as.matrix(met2)
# rownames(met2)=met2[,1]
# met2=met2[,2:ncol(met2)]
# dimnames3=list(rownames(met2),colnames(met2))
# met2=matrix(as.numeric(as.matrix(met2)),nrow=nrow(met2),dimnames=dimnames3)
# met2=avereps(met2)
# group=sapply(strsplit(colnames(met2),"\\-"), "[", 4)
# group=sapply(strsplit(group,""), "[", 1)
# group=gsub("2", "1", group)
# table(group)
# met2=met2[,group==0]
# rownames(met2)
# rownames(met) <- c("NID1", "NID2")
# rownames(met)
# met3 <- merge(met, met2, by= "row.names", all.x=TRUE)
# rownames(met3) <- met3[, 1]
# met3 <- met3[, -1]
# met <- met3
# 
# head(met)
# rownames(met)
# colnames(RNA)
# 
# methy <- met
# RNA <- data3
# rownames(methy)=paste(rownames(methy),"methy",sep="|")
# rownames(RNA)=paste(rownames(RNA),"exp",sep="|")
# sameSample=intersect(colnames(methy),colnames(RNA))
# sameSample
# length(sameSample)
# 
# merge=rbind(id=sameSample,methy[,sameSample],RNA[,sameSample])
# merge[1:4, 1:4]
# # write.table(merge, file="NID2_met_exp_LGG_merge.txt",sep="\t",quote=F,col.names=F)   
# write.table(merge, file="NID2_met_exp_GBM_merge.txt",sep="\t",quote=F,col.names=F) 


####
inputFile1 <- "../data/NID2_met_exp_LGG_merge.txt"
inputFile2 <- "../data/NID2_met_exp_GBM_merge.txt"
methyGene="NID2|methy"
expGene="NID2|exp"
rt1=read.table(inputFile1,sep="\t",header=T,check.names=F,row.names=1)
rt2=read.table(inputFile2,sep="\t",header=T,check.names=F,row.names=1)
sameSample=intersect(colnames(rt1),colnames(rt2))
length(sameSample)
rt1[1:4, 1:4]
type <- c(rep("LGG (N=529)", 529), rep("GBM (N=138)", 138))
exp <- cbind(rt1[expGene,], rt2[expGene,])
met <- cbind(rt1[methyGene,], rt2[methyGene,])
df1 <- data.frame(as.numeric(exp), log2(as.numeric(met)+1))
colnames(df1) <- c("NID2 Expression", "NID2 Methylation")
df1$Group <- type
table(df1$Group)
ordered.item <- c("LGG (N=529)", "GBM (N=138)")
ordered.item <- factor(1:length(ordered.item),labels = ordered.item)
df1$Group <-factor(df1$Group, levels = levels(ordered.item))
factor(df1$Group)

####
library(psych)
library("ggExtra")
library("ggpubr")
## main p
sp <- ggscatter(df1, x = "NID2 Expression", y = "NID2 Methylation",
          add = "reg.line",                         
          conf.int = TRUE,                         
          color = "Group", palette = "jco") +
  stat_cor(
    aes(label = paste(..r.label.., ..p.label.., "Gliomas", sep = "~`,`~")),
    label.x = 3 , show.legend = FALSE, cex = 5.5) 
sp
xplot <- ggdensity(df1, "NID2 Expression", fill = "Group",
                   palette = "jco")
yplot <- ggdensity(df1, "NID2 Methylation", fill = "Group", 
                   palette = "jco")+
  rotate()
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")
####
library(cowplot)
plot_grid(xplot, NULL, sp, yplot, ncol = 2, align = "hv", 
          rel_widths = c(2, 1), rel_heights = c(1, 2))

xdens <- axis_canvas(sp, axis = "x")+
  geom_density(data = df1, aes(x = `NID2 Expression`, fill = Group),
               alpha = 0.7, size = 0.2)+
  ggpubr::fill_palette("jco")
xdens
ydens <- axis_canvas(sp, axis = "y", coord_flip = TRUE)+
  geom_density(data = df1, aes(x = `NID2 Methylation`, fill = Group),
               alpha = 0.7, size = 0.2)+
  coord_flip()+
  ggpubr::fill_palette("jco")
p1 <- insert_xaxis_grob(sp, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")

pdf(file="Supplementary_Figure_2.pdf",width=10, height=10)
ggdraw(p2)
dev.off()

