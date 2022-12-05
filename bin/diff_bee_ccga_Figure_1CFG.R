
library(beeswarm)

## cgga
## survival_cgga_Figure_3EFGHIJ.R
head(clinic)
table(clinic$Grade)
rt <- clinic
table(rt$Histology)
rt$Histology <- gsub("rGBM", "GBM", rt$Histology)
rt$Histology <- gsub("sGBM", "GBM", rt$Histology)
rt$Histology <- gsub("rAOA", "LGG", rt$Histology)
rt$Histology <- gsub("AOA", "LGG", rt$Histology)
rt$Histology <- gsub("rAA", "LGG", rt$Histology)
rt$Histology <- gsub("rAO", "LGG", rt$Histology)
rt$Histology <- gsub("rOA", "LGG", rt$Histology)
rt$Histology <- gsub("AA", "LGG", rt$Histology)
rt$Histology <- gsub("AO", "LGG", rt$Histology)
rt$Histology <- gsub("OA", "LGG", rt$Histology)
rt$Histology <- gsub("rA", "LGG", rt$Histology)
rt$Histology <- gsub("rO", "LGG", rt$Histology)
rt$Histology <- gsub("O", "LGG", rt$Histology)
rt$Histology <- gsub("A", "LGG", rt$Histology)
rt$Histology
exp <- t(exp)
rownames(exp)
rownames(rt) <- rt$CGGA_ID
sample <- intersect(rownames(exp), rownames(rt))
df <- merge(exp[sample, ], rt[sample, ], by= "row.names")
colnames(df)

####
df <- as.data.frame(temp_ggplot)  #### GSE16011
df <- rt    #### GSE7696 and  GSE4290 start
head(df)
df <- df[, c("Row.names", "NID2", "Histology")] 
df <- df[, c("Row.names", "NID2", "group")]
colnames(df)
colnames(df)[3] <- "Type"
colnames(df)[2] <- "expression"
colnames(df) <- c("expression", "Type")
yMin=min(df$expression)                    
yMax=max(df$expression)
table(df$Type)
head(df)

wilcoxTest=wilcox.test(as.numeric(expression) ~ Type, data=df)
pvalue=wilcoxTest$p.value
pvalue
if (pvalue<0.001) { pvalue = "***" 
} else if (pvalue<0.01) { pvalue = "**" 
}else if (pvalue<0.05) { pvalue = "*" 
} else { pvalue=paste0("Pvalue = ",sprintf("%.03f",pvalue))
}
table(df$Type)
labels=c('LGG\n(N=625)','GBM\n(N=388)')     ## CGGA
labels=c('Non-tumor brain tissues\n(N=4)','GBM\n(N=80)')     ## GEO GSE7696
# labels=c('LGG\n(N=76)','GBM\n(N=77)')     ## GEO GSE4290
labels=c('Control\n(N=8)','Glioma\n(N=276)')     ## GEO GSE16011
labels=c('non-tumor\n(N=23)','GBM\n(N=77)')     ## GEO GSE4290
df$Type <- factor(df$Type, c("LGG", "GBM"))
df$Type <- factor(df$Type, c("control", "glioma"))
df$Type <- factor(df$Type, c("non-tumor", "GBM"))
df$expression <- as.numeric(df$expression)
pdf(file = "Figure_1F.pdf", width=6, height=5)
pdf(file = "Figure_1G.pdf", width=6, height=5)
pdf(file = "Figure_1C.pdf", width=6, height=5)
pdf(file = "Figure_3A.pdf", width=6, height=5)
par(mar=c(2,3,2,0.5), mgp=c(2, 1, 0), tck=0.01)
gene <- "NID2"
boxplot(expression ~ Type, data = df, col = "white", names=labels, xlab="",
        ylab = paste(gene," expression",sep=""), cex.axis=1, main = "")  ## final
stripchart(expression ~ Type,
           data = df,
           method = "jitter",
           pch = 19,
           col = c("#2b83ba","#d7191c"),  ## cgga
           # col = c("dark green","#d7191c"),
           vertical = TRUE,
           add = TRUE)
ySeg=yMax*1.01
segments(1,ySeg, 2,ySeg);segments(1,ySeg, 1,ySeg*0.98);segments(2,ySeg, 2,ySeg*0.98)
text(1.5,ySeg*1.005,labels=pvalue,cex=1.5)  ## final
dev.off()

