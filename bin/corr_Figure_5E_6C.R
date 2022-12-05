## survival_tcga_Figure_2EFGH.R
data[1:4, 1:4]
dim(data)
ins <- c("VEGFA", "FLT4", "KDR", "FLT1", "VWF", "PECAM1", "ANGPT1", "ANGPT2", "CDH5", "NID2")
exp <-data[row.names(data) %in% ins, ]
rownames(exp)


## angio
library(tidyr)
library(dplyr)
library(tibble)
rownames(exp)
exp[1:4, 1:4]
exp <- exp[sort(rownames(exp)), ]
rownames(exp)
rownames(exp) <- c("Angiopoietin (ANGPT) 1", "Angiopoietin (ANGPT) 2", "Vascular endothelial cadherin (CDH5)", "VEGFR1 (FLT1)", "VEGFR3 (FLT4)", "VEGFR2 (KDR)", "NID2", "CD31 (PECAM1)", "VEGFA", "von Willebrand factor (VWF)")
exp_temp <- as.data.frame(t(exp))
exp_longer<- exp_temp %>% 
  rownames_to_column(var = 'Sample') %>% 
  pivot_longer( cols =   c("Angiopoietin (ANGPT) 1", "Angiopoietin (ANGPT) 2", "Vascular endothelial cadherin (CDH5)", "VEGFR1 (FLT1)", "VEGFR3 (FLT4)", "VEGFR2 (KDR)", "CD31 (PECAM1)", "VEGFA", "von Willebrand factor (VWF)"),
                names_to = 'Gene',
                values_to = 'Expression')
head(exp_longer)

set.seed(123)

library(dplyr)
library(ggplot2)
library(ggpubr)
p <- ggplot(exp_longer, aes(x = NID2, y = Expression)) +
  geom_point() + 
  geom_smooth(method = "lm", col = "blue") +
  facet_wrap(~ Gene, nrow = 2) +
  stat_cor(method = "pearson", label.y = 9.5) +
  ylab("The expression of vascular-related marker genes") +
  xlab("NID2 expression")
pdf("Figure_5E.pdf",height=5, width=14)
print(p)
dev.off()

## immue
ins <- c("NID2", "PDCD1", "CD274", "IL12A", "IL12B", "CD28", "TGFB2", "TGFB1", "IL10", "CSF1", "STAT2", "PDCD1LG2", "CTLA4", "STAT3", "IL1B", "IL1A", "IL1R1", "IL1RAP", "IL1RN", "IL1R2", "IL2", "TNF", "IFNG", "IL6") 
exp <-data[row.names(data) %in% ins, ]
exp <- exp[sort(rownames(exp)), ]
rownames(exp)
rownames(exp) <- c("PD-L1 (CD274)", "CD28", "M-CSF (CSF1)", "CTLA4", "IFNG", "IL-10 (IL10)", "IL-12 (IL12A)", "IL-12 (IL12B)", "IL-1 (IL1A)", "IL-1 (IL1B)", "IL-1R (IL1R1)", "IL-1R (IL1R2)", "IL-1R (IL1RAP)", "IL-1R (IL1RN)", "IL-2 (IL2)", "IL-6 (IL6)","NID2",  "PD-1 (PDCD1)", "PD-L2 (PDCD1LG2)", "STAT2", "STAT3", "TGFB1", "TGFB2", "TNF")
exp[1:4, 1:4]
exp_temp <- as.data.frame(t(exp))
set.seed(123)
head(exp_temp)
exp_temp <- exp_temp %>% select(NID2, everything())

gene <- "NID2"
outTab=data.frame()
for(i in colnames(exp_temp)[2:(ncol(exp_temp))]){
  x=as.numeric(exp_temp[,gene])
  y=as.numeric(exp_temp[,i])
  if(sd(y)==0){y[1]=0.00001}
  cor=cor.test(x, y, method="spearman", exact=FALSE)
  outVector=cbind(Gene=i, cor=cor$estimate, pvalue=cor$p.value)
  outTab=rbind(outTab,outVector)
  if(cor$p.value<0.05){
    outFile=paste0("cor.", i, ".pdf")
    df1=as.data.frame(cbind(x,y))
    p1=ggplot(df1, aes(x, y)) + 
      xlab(paste0(gene, " expression")) + ylab(paste0(i, " expression"))+
      geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
      stat_cor(method = 'spearman', aes(x =x, y =y))
    pdf(file=outFile, width=4, height=3)
    print(p1)
    dev.off()
  }
}
outTab[outTab$Gene == "CTLA4", ]
write.table(outTab,file="cor.result.txt",sep="\t",row.names=F,quote=F)


x=as.numeric(exp_temp[, "NID2"])
y=as.numeric(exp_temp[, "STAT3"])
if(sd(y)==0){y[1]=0.00001}
cor=cor.test(x, y, method="spearman", exact=FALSE)
outVector=cbind(Gene=i, cor=cor$estimate, pvalue=cor$p.value)
outTab=rbind(outTab,outVector)
df1=as.data.frame(cbind(x,y))
p1=ggplot(df1, aes(x, y)) + 
  xlab(paste0(gene, " expression")) + ylab(paste0("STAT3", " expression"))+
  geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =x, y =y), label.y = 8)
pdf(file="cor.STAT3.pdf", width=4, height=3)
print(p1)
dev.off()

