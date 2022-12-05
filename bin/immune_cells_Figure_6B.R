library(limma)
library(reshape2)
library(ggpubr)
library(vioplot)
library(ggExtra)

immFile="CIBERSORT-Results.txt" 
pFilter=0.05

## survival_tcga_Figure_2EFGH.R
gene = "NID2"
data[1:4, 1:4]
data=avereps(data)
data=data[rowMeans(data)>0.5,]
data=t(data)
data=t(avereps(data))

geneExp=as.data.frame(t(data[gene,,drop=F]))
cutpoint <- "2.513645"   ## log2(tpm + 1)
# cutpoint <- "1.163174"
geneExp$Type=ifelse(geneExp[, gene] > cutpoint, "High", "Low")
head(geneExp)

immune=read.table(immFile, header=T, sep="\t", check.names=F)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])

rownames(immune) <- immune[,1]
immune=immune[,2:ncol(immune)]
dimnames1=list(rownames(immune),colnames(immune))
immune=matrix(as.numeric(as.matrix(immune)),nrow=nrow(immune),dimnames=dimnames1)
immune=avereps(immune)
immune[1:4, 1:4]

sameSample=intersect(row.names(immune), row.names(geneExp))
sameSample
rt=cbind(immune[sameSample,,drop=F], geneExp[sameSample,,drop=F])
rownames(rt)

data_imm=rt[,-(ncol(rt)-1)]
rownames(data_imm)
data_imm=melt(data_imm,id.vars=c("Type"))
colnames(data_imm)=c("Type", "Immune", "Expression")
head(data_imm)
group=levels(factor(data_imm$Type))
data_imm$gene=factor(data_imm$Type, levels=c("Low","High"))

bioCol <-  c("#2b83ba","#d7191c")
boxplot=ggboxplot(data_imm, x="Immune", y="Expression", fill="gene",
                  xlab="",
                  ylab="Fraction",
                  legend.title=gene,
                  width=0.8,
                  palette=bioCol)+
  rotate_x_text(50)+
  stat_compare_means(aes(group=gene),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")

pdf(file="Figure_6B.pdf", width=9, height=6)
print(boxplot)
dev.off()


outTab=data.frame()
for(i in colnames(rt)[1:(ncol(rt)-2)]){
  x=as.numeric(rt[,gene])
  y=as.numeric(rt[,i])
  if(sd(y)==0){y[1]=0.00001}
  cor=cor.test(x, y, method="spearman")
  outVector=cbind(Cell=i, cor=cor$estimate, pvalue=cor$p.value)
  outTab=rbind(outTab,outVector)
  if(cor$p.value<0.05){
    outFile=paste0("cor.", i, ".pdf")
    df1=as.data.frame(cbind(x,y))
    p1=ggplot(df1, aes(x, y)) + 
      xlab(paste0(gene, " expression")) + ylab(i)+
      geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
      stat_cor(method = 'spearman', aes(x =x, y =y))
    # pdf(file=outFile, width=5.2, height=5)
    # print(p1)
    # dev.off()
  }
}

write.table(outTab,file="cor.result.txt",sep="\t",row.names=F,quote=F)
