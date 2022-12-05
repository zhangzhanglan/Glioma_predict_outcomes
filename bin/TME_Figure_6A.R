
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")


library(limma)
library(reshape2)
library(ggpubr)

scoreFile="estimateScores.txt"
## survival_tcga_Figure_2EFGH.R 

gene <- "NID2"
temp <- t(exp)
temp=avereps(temp)
exp <- temp
# exp$Type=ifelse(exp[,gene]>median(exp[,gene]), "High", "Low")
# cutpoint <- "1.163174"
cutpoint <- "2.513645"   ## log2(tpm + 1)
exp <- as.data.frame(exp)
exp[, gene]
exp$Type=ifelse(exp[, gene] > cutpoint, "High", "Low")
head(exp)

score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
score=score[,1:3]
colnames(score) <- c("Stromal Score", "Immune Score", "ESTIMATE Score")

sameSample=intersect(row.names(exp), row.names(score))
sameSample
exp=exp[sameSample,"Type",drop=F]
score=score[sameSample,,drop=F]
rt=cbind(score, exp)
rt$Type=factor(rt$Type, levels=c("Low", "High"))

data=melt(rt, id.vars=c("Type"))
colnames(data)=c("Type", "scoreType", "Score")


p=ggviolin(data, x="scoreType", y="Score", fill = "Type",
	     xlab="",
	     ylab="TME score",
	     legend.title=gene,
	     add = "boxplot", add.params = list(color="white"),
	     palette = c("#2b83ba","#d7191c"), width=1)
# p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Type),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif") + rotate_x_text(50)

pdf(file="Figure_6A.pdf", width=5, height=5)
print(p1)
dev.off()
