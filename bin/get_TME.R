library(limma)
library(estimate)

# survival_tcga_Figure_2EFGH.R
outTab=data.frame()
data[1:4, 1:4]
out=data[rowMeans(data)>0,]
out=rbind(ID=colnames(out),out)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)
	
filterCommonGenes(input.f="uniq.symbol.txt", output.f="commonGenes.gct", id="GeneSymbol")
estimateScore(input.ds = "commonGenes.gct", output.ds="estimateScore.gct")
	
scores=read.table("estimateScore.gct",skip = 2,header = T)
scores[1:4, 1:4]
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)
rownames(scores)=gsub("\\.","\\-",rownames(scores))
file.remove("commonGenes.gct")
file.remove("estimateScore.gct")
# file.remove("uniq.symbol.txt")

out=cbind(ID=row.names(scores),scores)
write.table(out,file="estimateScores.txt",sep="\t",quote=F,row.names=F)

