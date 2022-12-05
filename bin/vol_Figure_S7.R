library("limma")

rt2=read.table("tcgaSymbol_LGG.txt",sep="\t",header=T,check.names=F)
rt2=as.matrix(rt2)
rownames(rt2)=rt2[,1]
exp2=rt2[,2:ncol(rt2)]
dimnames2=list(rownames(exp2),colnames(exp2))
data2=matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)
data2=avereps(data2)
data2=data2[rowMeans(data2)>0.5,]
dim(data2)
group=sapply(strsplit(colnames(data2),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
table(group)
data2=data2[,group==0]
dim(data2)

rt3=read.table("tcgaSymbol_GBM.txt",sep="\t",header=T,check.names=F)
rt3=as.matrix(rt3)
rownames(rt3)=rt3[,1]
exp3=rt3[,2:ncol(rt3)]
dimnames3=list(rownames(exp3),colnames(exp3))
data3=matrix(as.numeric(as.matrix(exp3)),nrow=nrow(exp3),dimnames=dimnames3)
data3=avereps(data3)
data3=data3[rowMeans(data3)>0.5,]
dim(data3)
group=sapply(strsplit(colnames(data3),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
table(group)
data3=data3[,group==0]
dim(data3)

dim(data2)
dim(data3)
sameGene=intersect( row.names(data2), row.names(data3))
data=cbind(data2[sameGene,], data3[sameGene,])

ccga=read.table("all_genes_ccga.txt", sep="\t",header=T,check.names=F)
data <- data[which(row.names(data) %in% ccga$x), ]
dim(data)

rownames(data)
ins <- data[rownames(data) == "NID2", ]
ins <- as.data.frame(ins)
colnames(ins) <- "gene"
head(ins)
median(ins$gene)
sur.cut$cutpoint$cutpoint
ins$grade <- ifelse(ins$gene <= 1.163174,1,2)
head(ins)
ins <- ins[order(ins$grade),]
table(ins$grade)

dataL=data[, row.names(ins[ins[,"grade"]==1,])]
dataH=data[, row.names(ins[ins[,"grade"]==2,])]
dim(dataL)
dim(dataH)
dim(datar)
datar=cbind(dataL, dataH)
head(datar)
table(ins$grade)

conNum=441
treatNum=256

outTab=data.frame()
grade=c(rep(1,conNum),rep(2,treatNum))
dim(datar)
table(grade)
for(i in row.names(datar)){
  geneName=unlist(strsplit(i,"\\|",))[1]
  geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=datar[i,], grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
  conGeneMeans=mean(datar[i,1:conNum])
  treatGeneMeans=mean(datar[i,(conNum+1):ncol(datar)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  pvalue
  logFC
  conMed=median(datar[i,1:conNum])
  treatMed=median(datar[i,(conNum+1):ncol(datar)])
  diffMed=treatMed-conMed
  diffMed
	# if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){
		  outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
	 # }
}

fdrFilter=0.001
logFCfilter=1

pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)
head(outTab)
write.table(outTab,file="all_bestcutoff.txt",sep="\t",row.names=F,quote=F)
outDiff <- subset(outTab, (abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter))
outDiff[outDiff$gene == "NID2", ]
write.table(outDiff,file="diff_bestcutoff.xls",sep="\t",row.names=F,quote=F)
dim(outDiff)

heatmap=rbind(ID=colnames(datar[as.vector(outDiff[,1]),]),datar[as.vector(outDiff[,1]),])
write.table(heatmap,file="diffGeneExp_bestcutoff.txt",sep="\t",col.names=F,quote=F)

# pdf(file="vol_bestcutoff.pdf",height=5,width=5)
# xMax=max(abs(as.numeric(as.vector(outTab$logFC))))
# yMax=max(-log10(outTab$fdr))+1
# plot(as.numeric(as.vector(outTab$logFC)), -log10(outTab$fdr), xlab="logFC",ylab="-log10(fdr)",
#      main="DEGs of NID2 Expression (High VS Low)", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=0.8)
# diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))>logFCfilter)
# points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="#d7191c",cex=0.8)
# diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))<(-logFCfilter))
# dim(diffSub)
# points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="#fdae61",cex=0.8)
# legend("topleft", c("Down", "Not", "Up"), fill = c('#fdae61', 'black', '#d7191c'), title= "Group", cex=0.6)
# text(-1, 60, "FDR adjusted p < 0.001", cex=0.6)
# text(-1, 55, "|logFC| > 1", cex=0.6)
# abline(v=0,lty=2,lwd=3)
# dev.off()

## final
library(ggplot2)
library(dplyr)
pdf(file="vol_bestcutoff.pdf",height=5,width=14)
head(data_v)
data_v <- outTab
data_v$logFC <- as.numeric(as.vector(data_v$logFC))
data_v$fdr <- as.numeric(as.vector(data_v$fdr))

data_v$Group = ifelse(data_v$fdr < fdrFilter & abs(data_v$logFC) >= logFCfilter, 
                        ifelse(data_v$logFC> logFCfilter ,'Up','Down'),
                        'Stable')
table(data_v$Group)
for_label <- data_v %>%
  filter(abs(logFC) >1.5 & -log10(fdr) > -log10(0.000001))

p <- ggplot(data = data_v, 
            aes(x = logFC, 
                y = -log10(fdr))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=Group)) +
  scale_color_manual(values=c("#2b83ba", "grey","#d7191c"), guide = guide_legend(override.aes = list(size = 8)))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  coord_flip() +
  theme_bw() + theme(panel.grid=element_blank()) 

p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = gene),
    data = for_label,
    color="black",
    box.padding = 0.5, max.overlaps = Inf
  )

dev.off()

library(pheatmap)
hmExp=data[as.vector(outDiff[,1]),]
hmExp=log2(hmExp+0.001)
Type=c(rep("NID2 Expression Low",conNum),rep("NID2 Expression High",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)


pdf(file="heatmap.pdf",height=12,width=15)
pheatmap(hmExp, 
         annotation=Type, 
         # color = colorRampPalette(c("#fdae61", "black", "#d7191c"))(10),
         # color = heat.colors(8),
         cluster_cols =F,
         cluster_row =F,
         show_colnames = F,
         show_rownames = F,
         fontsize = 12,
         fontsize_row=3,
         fontsize_col=10)
dev.off()
  
