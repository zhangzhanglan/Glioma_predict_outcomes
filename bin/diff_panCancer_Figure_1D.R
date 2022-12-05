library(ggpubr)

data=read.table("diff_panCancer_singleGeneExp.txt", sep="\t", header=T, check.names=F)
gene=colnames(data)[2]
colnames(data)[2]="expression"
table(data$CancerType)
data$CancerType <- factor(data$CancerType, levels = c("GBM", "LGG", "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "HNSC", "KICH",  "KIRC", "KIRP", "LAML", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA",  "THYM",  "UCEC", "UCS", "UVM"))
p=ggboxplot(data, x="CancerType", y="expression", color = "Type",
     ylab=paste0(gene," expression"),
     xlab="",
     palette = c("dark green","#d7191c") )
p=p+rotate_x_text(60)
pdf(file="Figure_1D.pdf", width=8,height=5)
p+stat_compare_means(aes(group=Type),
      method="wilcox.test",
      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
      label = "p.signif")
dev.off()
factor(levels(data$CancerType))
# data[data$CancerType == "LGG", ]
outTab=data.frame()
for (i in factor(levels(data$CancerType))) {
      tab1=table(data[data$CancerType == i, "Type"])
      labelNum=length(tab1)
      if(labelNum==1){next }
      res = compare_means(expression ~ Type, data = data[data$CancerType == i, ])
      outTab=rbind(outTab, cbind(i, res))
}
outTab
