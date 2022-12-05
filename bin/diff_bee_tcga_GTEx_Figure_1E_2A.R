library("limma")

inputFile="merge_NID2.txt"
conNum = 1152+5
treatNum = 168+529

outTab=data.frame()
grade=c(rep(1, conNum), rep(2, treatNum))
table(grade)
rt=read.table(inputFile, sep="\t", header=T,check.names=F)
rt=as.matrix(rt)
head(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
exp[1:2, 1:2]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

library(beeswarm)
gene="NID2"

normal=1157
GBM=168
LGG=529
Type=c(rep(1, normal), rep(2, GBM), rep(3, LGG))
nor_exp <- data[, 1:1157]
gbm_exp <- data[, 1158:1325]
lgg_exp <- data[, 1326:1854]
data_nor_gbm <- cbind(nor_exp, gbm_exp)
data_nor_lgg <- cbind(nor_exp, lgg_exp)
data_lgg_gbm <- cbind(lgg_exp, gbm_exp)

type_nor_gbm <- c(rep(1,normal), rep(2,GBM))
type_nor_lgg <- c(rep(1,normal), rep(2,LGG))
type_lgg_gbm <-  c(rep(2,LGG), rep(1,GBM))

## nor_lgg
# pdf(file="diff_tcga_normal_VS_LGG.pdf", width=6, height=5)
par(mar=c(2,3,2,0.5), mgp=c(2, 1, 0), tck=0.01)
rt=rbind(expression=data_nor_lgg[gene,],Type=type_nor_lgg)
rt=as.matrix(t(rt))
wilcoxTest=wilcox.test(expression ~ Type, data=rt)
pvalue=wilcoxTest$p.value
pvalue
if(pvalue<0.001){ pvalue = "***"
}else{ pvalue=paste0("=" , sprintf("%.03f",pvalue))}
pvalue
df <- as.data.frame(rt)
yMin=min(df$expression)                    
yMax=max(df$expression)
table(df$Type)
df$Type <- factor(df$Type, c("1", "2"))
labels=c('Non-tumor brain tissues\n(N=1157)','LGG\n(N=529)')
boxplot(expression ~ Type, data = df, col = "white", names=labels, xlab="",
        ylab = paste(gene," expression",sep=""), cex.axis=1)
stripchart(expression ~ Type,
           data = df,
           method = "jitter",
           pch = 19,
           col = c("dark green","#2b83ba"),
           vertical = TRUE,
           add = TRUE)
ySeg=yMax*1.01
segments(1,ySeg, 2,ySeg);segments(1,ySeg, 1,ySeg*0.98);segments(2,ySeg, 2,ySeg*0.98)
text(1.5,ySeg*1.005,labels=pvalue,cex=1.5)
# dev.off()
## nor_gbm
pdf(file="Figure_1E.pdf", width=6, height=5)
par(mar=c(2,3,2,0.5), mgp=c(2, 1, 0), tck=0.01)
rt=rbind(expression=data_nor_gbm[gene,],Type=type_nor_gbm)
rt=as.matrix(t(rt))
wilcoxTest=wilcox.test(expression ~ Type, data=rt)
pvalue=wilcoxTest$p.value
pvalue
if( pvalue<0.001){ pvalue = "***"
}else{ pvalue=paste0("=",sprintf("%.03f",pvalue)) }
df <- as.data.frame(rt)
yMin=min(df$expression)                    
yMax=max(df$expression)
table(df$Type)
df$Type <- factor(df$Type, c("1", "2"))
labels=c('Non-tumor brain tissues\n(N=1157)','GBM\n(N=168)')
boxplot(expression ~ Type, data = df, col = "white", names=labels, xlab="", ylab = paste(gene," expression",sep=""), cex.axis=1)
stripchart(expression ~ Type,
           data = df,
           method = "jitter",
           pch = 19,
           col = c("dark green","#d7191c"),
           vertical = TRUE,
           add = TRUE)
ySeg=yMax*1.01
segments(1,ySeg, 2,ySeg);segments(1,ySeg, 1,ySeg*0.98);segments(2,ySeg, 2,ySeg*0.98)
text(1.5,ySeg*1.005,labels=pvalue,cex=1.5)
dev.off()
## lgg_gbm
pdf(file="Figure_2A.pdf", width=6, height=5)
par(mar=c(2,3,0.5,0.5), mgp=c(2, 1, 0), tck=0.01)
rt=rbind(expression=data_lgg_gbm[gene,],Type=type_lgg_gbm)
rt=as.matrix(t(rt))
wilcoxTest=wilcox.test(expression ~ Type, data=rt)
pvalue=wilcoxTest$p.value
pvalue
if( pvalue<0.001){ pvalue = "***"
}else{ pvalue=paste0("=",sprintf("%.03f",pvalue)) }
df <- as.data.frame(rt)
yMin=min(df$expression)                    
yMax=max(df$expression)
table(df$Type)
df$Type <- factor(df$Type, c("2", "1"))
labels=c('LGG\n(N=529)','GBM\n(N=168)')
boxplot(expression ~ Type, data = df, col = "white", names=labels, xlab="", ylab = paste(gene," expression",sep=""), cex.axis=1)
stripchart(expression ~ Type,
           data = df,
           method = "jitter",
           pch = 19,
           col = c("#2b83ba","#d7191c"),
           vertical = TRUE,
           add = TRUE)
ySeg=yMax*1.01
segments(1,ySeg, 2,ySeg);segments(1,ySeg, 1,ySeg*0.98);segments(2,ySeg, 2,ySeg*0.98)
text(1.5,ySeg*1.005,labels=pvalue,cex=1.5)
dev.off()

####
df <- t(data)
Type=c(rep("Normal",normal), rep("GBM",GBM), rep("LGG", LGG))
df <- as.data.frame(df)
df$Type <- Type
head(df)
library(ggpubr)  
dim(data)
p=ggboxplot(df, x="Type", y="NID2", color = "Type",
            ylab=paste0(gene," expression"),
            xlab="",
            palette = c("#abdda4","#d7191c", "#2b83ba") )
p+stat_compare_means(aes(group=Type),
          # method="wilcox.test",
          method = "anova", label.y = 7
          # symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
          # label = "p.signif")
)

head(rt)
df <- as.data.frame(t(rt))
colnames(df) <- df[1, ]
df <- df[-1, ]
head(df)
df$Type <- grade
head(df)
df <- df[, -1]
colnames(df)[1] <- "expression"
par(mar=c(2,3,0.5,0.5), mgp=c(2, 1, 0), tck=0.01)
head(df)
rt=as.matrix(df)
wilcoxTest=wilcox.test(as.numeric(expression) ~ Type, data=rt)
pvalue=wilcoxTest$p.value
pvalue
if(pvalue<0.001){ pvalue = "***"
}else{ pvalue=paste0("=" , sprintf("%.03f",pvalue))}
table(grade)
gene <- "NID2"
labels=c('Non-tumor brain tissues\n(N=1157)','Gliomas\n(N=697)')
df$expression <- as.numeric(df$expression)
boxplot(expression ~ Type, data = df, col = "white", names=labels, xlab="",
        ylab = paste(gene," expression",sep=""), cex.axis=1)
stripchart(expression ~ Type,
           data = df,
           method = "jitter",
           pch = 19,
           col = c("dark green","#d7191c"),
           vertical = TRUE,
           add = TRUE)
yMin=min(df$expression)                    
yMax=max(df$expression)
ySeg=yMax*1.01
segments(1,ySeg, 2,ySeg);segments(1,ySeg, 1,ySeg*0.98);segments(2,ySeg, 2,ySeg*0.98)
text(1.5,ySeg*1.005,labels=pvalue,cex=1.5)

