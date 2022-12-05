library(limma)

rt1=read.table("../data/GTExNormalBrainExp.txt",sep="\t",header=T,check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
data1=matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames1)
data1=avereps(data1)

rt2=read.table("../data/tcgaSymbol_LGG.txt",sep="\t",header=T,check.names=F)
rt2=as.matrix(rt2)
rownames(rt2)=rt2[,1]
exp2=rt2[,2:ncol(rt2)]
dimnames2=list(rownames(exp2),colnames(exp2))
data2=matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)
data2=avereps(data2)

rt3=read.table("../data/tcgaSymbol_GBM.txt",sep="\t",header=T,check.names=F)
rt3=as.matrix(rt3)
rownames(rt3)=rt3[,1]
exp3=rt3[,2:ncol(rt3)]
dimnames3=list(rownames(exp3),colnames(exp3))
data3=matrix(as.numeric(as.matrix(exp3)),nrow=nrow(exp3),dimnames=dimnames3)
data3=avereps(data3)

dim(data1)
head(colnames(data1))

sameGene1=intersect( row.names(data1), row.names(data3))
sameGene <- intersect(sameGene1,  row.names(data2))
data=cbind(data1[sameGene,], data3[sameGene,], data2[sameGene,])

conNum=1152+5
treatNum = 168+529
grade=c(rep(1,conNum),rep(2,treatNum))
design <- as.factor(grade)
design <- model.matrix(~0 + design)
batch <- c(rep(1,1152),rep(2, 702))
table(batch)
data_rm <- removeBatchEffect(data, batch = batch, design = design)
outTab <- data_rm

outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab,file="merge_rb.txt",sep="\t",quote=F,col.names=F)

## head -1 merge_rb.txt > merge_NID2.txt
## grep -w NID1 merge_rb.txt >>  merge_NID2.txt
## grep -w NID2 merge_rb.txt >>  merge_NID2.txt
