library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(corrplot)

## survival_tcga_Figure_2EFH.R
pFilter=0.001   
geneName="NID2"
geneFile="../data/CheckpointGene.txt"

data[1:4, 1:4]
data=avereps(data)

gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data), as.vector(gene[,1]))
sameGene
data=t(data[c(geneName, sameGene),])
data=t(avereps(data))

x=as.numeric(data[geneName,])
outTab=data.frame()
for(i in sameGene){
	if(i==geneName){next}
    y=as.numeric(data[i,])
	corT=cor.test(x, y, method = 'pearson')
	cor=corT$estimate
	pvalue=corT$p.value
	if(pvalue<pFilter){
		outTab=rbind(outTab, cbind(Query=geneName, Gene=i, cor, pvalue))
	}
}
write.table(file="../data/corResult.txt", outTab, sep="\t", quote=F, row.names=F)

data=t(data[c(geneName, as.vector(outTab[,2])),])
M=cor(data)

col1=colorRampPalette(colors =c("red","white","darkgreen"), space="Lab") 
pdf("Supplementary_Figure_8B.pdf",height=6,width=8)
corrplot(M, type = "lower", method = "color", tl.cex=0.6, col = col1(10),
         p.mat = res1$p, sig.level = c(.001, .01, .05), outline="#A6CEE3",
         insig = "label_sig", pch.cex = 0.5, pch.col = "black")
dev.off()
#
