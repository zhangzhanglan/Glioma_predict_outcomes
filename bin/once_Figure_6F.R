rm(list = ls())

CVCl <- read.csv("/Users/lanzhzh/Downloads/expasy_brain_celllines.csv", row.names = 1, header = FALSE)
head(CVCl)
CVCl

library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
dir='/Users/lanzhzh/Downloads/DataFiles/DataFiles/Training\ Data'
CTRP2 <- readRDS(file=file.path(dir,'CTRP2_Expr (TPM, not log transformed).rds'))
CTRP2 <- log10(CTRP2+1)
GDSC2_Res <- readRDS(file = file.path(dir,"CTRP2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res)
table(rownames(GDSC2_Res))
list <- colnames(GDSC2_Res)
write.table(list, file="drug_w.txt",sep="\t",quote=F,col.names=F)  
# CVCl <- ins_cells
# CVCl_brain <- intersect(rownames(GDSC2_Res), rownames(CVCl))
# CVCl_brain <- intersect(rownames(GDSC2_Res), CVCl)    ## CVCl == 21
# GDSC2_Res <- GDSC2_Res[CVCl_brain, ]
# GDSC2_Res <- GDSC2_Res[CVCl_brain, ins]
GDSC2_Res <- GDSC2_Res[, ins]
# CTRP2 <- CTRP2[, CVCl_brain]
head(GDSC2_Res)
GDSC2_Res <- na.omit(GDSC2_Res)
CTRP2 <- CTRP2[, rownames(GDSC2_Res)]
removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
GDSC2_Res <- removeColsAllNa(GDSC2_Res)
table(colnames(CTRP2) %in% "CVCL_4773")
dim(CTRP2)
dim(GDSC2_Res)
ins_cells <- rownames(GDSC2_Res)
table(colnames(GDSC2_Res) %in% "ZM")
colnames(GDSC2_Res)


f_rm_duplicated <- function(NameL, reverse=F) {
  tmp <- data.frame(table(NameL))
  if (reverse) {
    tmp <- tmp$NameL[tmp$Freq > 1]
  }else {
    tmp <- tmp$NameL[tmp$Freq == 1]
  }
  which(NameL %in% as.character(tmp))
}

colnames(data)
dim(data)
tpm <- data
tpm[1:4, 1:4]
table(duplicated(colnames(tpm)) )

comm <- intersect(rownames(CTRP2), rownames(tpm))
CTRP2 <- CTRP2[comm,]
tpm <- tpm[comm,]
dim(tpm)
tpm[1:4, 1:4]
as.matrix(tpm)[1:4, 1:4]
library(oncoPredict)
keep <- rowSums(CTRP2) > 0.8*ncol(CTRP2)
dim(CTRP2)
dim(GDSC2_Res)
GDSC2_Res[, "CIL55"]

## start
testPtype <- read.csv('./calcPhenotype_Output/DrugPredictions.csv', row.names = 1)
testPtype <- log(testPtype)
testPtype
colnames(testPtype)
table(colnames(testPtype) %in% "bev-a-CIZ-oo-mab")
ins <- c("pazopanib", 
         "sunitinib",
         "sorafenib", 
         "axitinib",
         "etoposide",
         "dasatinib", 
         "vandetanib",
         "tivozanib",
         "cabozantinib",
         "lenvatinib",
         "regorafenib",
         "foretinib",
         "vandetanib",
         "brivanib",
         "cediranib",
         "nintedanib",
         "linifanib",
         "RAF265",
         "tivozanib",
         "brivanib",
         "temozolomide",
         "procarbazine")
ins_testPtype <- testPtype[, colnames(testPtype) %in% ins]
ins_testPtype[1:4, 1:4]

gene="NID2"

geneExp=as.data.frame(t(data[gene,,drop=F]))
# geneExp$Type=ifelse(geneExp[,gene]>median(geneExp[,gene]), "High", "Low")
cutpoint <- "0.7566827"  ## log10(tpm + 1)

geneExp$Type=ifelse(geneExp[, gene] > cutpoint, "High", "Low")
head(geneExp)

comm <- intersect(rownames(ins_testPtype), rownames(geneExp))

## turn to person
df <- cbind(ins_testPtype[comm, ], geneExp[comm, ]$Type)
colnames(df)[[ncol(df)]]  <- "type"
head(df)

outTab=data.frame()
library(ggpubr)
library(Hmisc)
seq(ncol(df)-1)
for (i in  seq(ncol(df)-1)){
  rt <- df[, c(i, ncol(df))]
  colnames(rt)[1]
  drug_name <- colnames(rt)[1]
  colnames(rt) <- c("drug", "type")
  test=wilcox.test(drug~type, data=rt)
  lowMeans = median(rt[rt$type == "Low", ]$drug)
  higMeans = median(rt[rt$type == "High", ]$drug)
  logFC=log2(higMeans)-log2(lowMeans)
  diffPvalue=test$p.value
  drug_name <- capitalize(drug_name)
  outTab=rbind(outTab,cbind(drug=drug_name, lowMeans=lowMeans, higMeans=higMeans,logFC=logFC,pValue=diffPvalue))
  my_comparisons <- list(c("Low", "High"))
  head(rt)
  if(diffPvalue < 0.05){
    h <- ggviolin(rt, x="type", y="drug", fill = "type", 
         legend.title="NID2",
         xlab="",
         ylab=paste0(drug_name, " senstivity (IC50)"),
         palette=c("#2b83ba","#d7191c"),
         add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons, symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
    pdf(file=paste0("durgSenstivity2.", drug_name, ".pdf"), width=5, height=4.5)
    print(h)
    dev.off()
  }
}
outTab

comm
head(geneExp)
colnames(geneExp)
rt <- cbind(testPtype[comm, ], geneExp[comm, ])
colnames(rt)[[ncol(rt)]]  <- "NID2"
head(rt)
outTab=data.frame()
gene <- "NID2"
colnames(rt)
for(i in colnames(rt)[1:(ncol(rt)-1)]){
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
  }
}
# outTab 
write.table(outTab,file="drug.result.txt",sep="\t",row.names=F,quote=F)
