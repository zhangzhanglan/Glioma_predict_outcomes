library(limma)
library(ggpubr)
library(sva)
library(limma)
file1="../data/tcgaSymbol_LGG.txt"
file2="../data/tcgaSymbol_GBM.txt"

rt1=read.table(file1,sep="\t",header=T,check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
data1=matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames1)
data1=avereps(data1)
group=sapply(strsplit(colnames(data1),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
table(group)
data1=data1[,group==0]
colnames(data1) <- gsub("-0[12][ABC]$", "", colnames(data1))
data1[1:4, 1:4]

rt2=read.table(file2,sep="\t",header=T,check.names=F)
rt2=as.matrix(rt2)
rownames(rt2)=rt2[,1]
exp2=rt2[,2:ncol(rt2)]
dimnames2=list(rownames(exp2),colnames(exp2))
data2=matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)
data2=avereps(data2)
group=sapply(strsplit(colnames(data2),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
table(group)
colnames(data2[,group==1])
data2=data2[,group==0]
colnames(data2) <- gsub("-0[12][ABC]$", "", colnames(data2))
# colnames(data2) <- gsub("-11A$", "", colnames(data2))
# colnames(data2) <- gsub("[ABC]$", "", colnames(data2))
data2[1:4, 1:4]

dim(data1)
dim(data2)

sameGene=intersect(row.names(data1),row.names(data2))
data=cbind(data1[sameGene,],data2[sameGene,])

dim(data)

table(duplicated(colnames(data)))
data <- t(data)
data <- avereps(data)
data <- t(data)
dim(data) 

tciaFile="../data/TCIA-ClinicalData.tsv"
head(geneExp)

data[1:4, 1:4]
data=avereps(data)
data=data[rowMeans(data)>0.5,]
data=t(data)
data=t(avereps(data))

gene = "NID2"
geneExp=as.data.frame(t(data[gene,,drop=F]))
# geneExp$Type=ifelse(geneExp[,gene]>median(geneExp[,gene]), "High", "Low")
cutpoint <- "1.163174"
geneExp$Type=ifelse(geneExp[, gene] > cutpoint, "High", "Low")
head(geneExp)

ips=read.table(tciaFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(ips)
table(ips$clinical_data_tumor_tissue_site)
ips <- ips[, c(15, 16, 17, 18)]
head(ips)
ips <- na.omit(ips)
sameSample=intersect(row.names(ips), row.names(geneExp))
sameSample
ips=ips[sameSample, , drop=F]
geneExp=geneExp[sameSample, "Type", drop=F]
data_tica=cbind(ips, geneExp)
head(data_tica)

group <- levels(factor(data_tica$Type))
group
data_tica$Type=factor(data_tica$Type, levels=c("Low", "High"))
group=levels(factor(data_tica$Type))
group
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
colnames(data_tica)[1:(ncol(data_tica)-1)]

data_tica[data_tica$Type =="Low", "ips_ctla4_neg_pd1_neg"]
mean(data_tica[data_tica$Type =="High", "ips_ctla4_neg_pd1_neg"])
mean(data_tica[data_tica$Type =="Low", "ips_ctla4_neg_pd1_neg"])

mean(data_tica[data_tica$Type =="High", "ips_ctla4_neg_pd1_pos"])
mean(data_tica[data_tica$Type =="Low", "ips_ctla4_neg_pd1_pos"])

mean(data_tica[data_tica$Type =="High", "ips_ctla4_pos_pd1_neg"])
mean(data_tica[data_tica$Type =="Low", "ips_ctla4_pos_pd1_neg"])

mean(data_tica[data_tica$Type =="High", "ips_ctla4_pos_pd1_pos"])
mean(data_tica[data_tica$Type =="Low", "ips_ctla4_pos_pd1_pos"])

library(tidyr)
library(dplyr)
library(tibble)
## final not final, still final
plot_names <- c('ctla4_neg_pd1_neg' = "Control",
                'ctla4_neg_pd1_pos' = "PD-1 blocker",
                'ctla4_pos_pd1_neg' = "CTLA-4 blocker",
                'ctla4_pos_pd1_pos' = "CTLA-4 and PD-1 blockers")
my_comparisons <- list( c("Low", "High") )
p <- ggviolin(data_tica_longer, x = "NID2", y = "Expression", 
                ylab="Immunophenoscore",
                color = "NID2", palette = c("#2b83ba","#d7191c"),
                add = "jitter") +
  facet_wrap(vars(Treatment), labeller = as_labeller(plot_names), nrow = 1)
p <- p + stat_compare_means(comparisons = my_comparisons,
                    method="wilcox.test",
                    symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")), label = "p.signif",
                    label.y = 10, bracket.size = 0.5, 
                    )

pdf(file="Figure_6E.pdf", width=4.8*3, height=4.25)
print(p)
dev.off()
