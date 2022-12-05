# install.packages('e1071')

# if (!requireNamespace("BiocManager", quietly = TRUE))
   # install.packages("BiocManager")
# BiocManager::install("preprocessCore")

library("limma")
# survival_tcga_Figure_2EFGH.R
data[1:4, 1:4]
data=avereps(data)
data=data[rowMeans(data)>0,]

v <-voom(data, plot = F, save.plot = F)
out=v$E
out=rbind(ID=colnames(out),out)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)       

source("get_immune_cells_run.R")
results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=100, QN=TRUE)

#
