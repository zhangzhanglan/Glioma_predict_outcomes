library("GSVA")
library(fgsea)
gmtFile="../data/immue_suppression_sets.txt"
pathways.hallmark <- gmtPathways(gmtFile)

pathways.hallmark$`Markers of Tregs`
sameGene <- intersect(pathways.hallmark$`Markers of Tregs`, rownames(data))
sameGene

pathways.hallmark$`Immunosuppressive signaling pathways`
sameGene <- intersect(pathways.hallmark$`Immunosuppressive signaling pathways`, rownames(data))
sameGene

pathways.hallmark$Immunosuppressors
sameGene <- intersect(pathways.hallmark$Immunosuppressors, rownames(data))
sameGene

pathways.hallmark$`Immunosuppressive cytokines and checkpoints`
sameGene <- intersect(pathways.hallmark$`Immunosuppressive cytokines and checkpoints`, rownames(data))
sameGene

pathways.hallmark$`Tumor-supportive macrophage chemotactic and skewing molecules`
sameGene <- intersect(pathways.hallmark$`Tumor-supportive macrophage chemotactic and skewing molecules`, rownames(data))
sameGene

setdiff(pathways.hallmark$`Immune activators`, sameGene)

# survival_tcga_Figure_2EFGH.R
data[1:4, 1:4]
data=avereps(data)
# data=data[rowMeans(data)>0.5,]   ## 
# data=t(data)
# data=t(avereps(data))
rt <- data
rownames(rt)
# leukemia_es <- gsva(rt, pathways.hallmark, min.sz=10, max.sz=500)
leukemia_es <- gsva(rt, pathways.hallmark, min.sz=4, max.sz=500)
rownames(leukemia_es)

gene="NID2"
geneExp=as.data.frame(t(data[gene,,drop=F]))
# geneExp$Type=ifelse(geneExp[,gene]>median(geneExp[,gene]), "High", "Low")
cutpoint <- "2.513645"   ## log2(tpm + 1)
# cutpoint <- "1.163174"
geneExp$Type=ifelse(geneExp[, gene] > cutpoint, "High", "Low")
head(geneExp)
geneExp <- geneExp[order(geneExp[, "NID2"]),]
rownames(geneExp)
colnames(leukemia_es)

###
sample <- intersect(rownames(geneExp), colnames(leukemia_es))
sample
geneExp <- geneExp[sample, ]
leukemia_es <- leukemia_es[, sample]
leukemia_es[1:4, 1:4]
colnames(geneExp)
colnames(geneExp)[2] <- "Group"
data["STAT5B", ]

library(pheatmap)
annotation_col = as.data.frame(geneExp[, -1])
rownames(annotation_col) <- rownames(geneExp)
colnames(annotation_col) <- "NID2"
head(annotation_col)
annotation_col$NID2 <- as.factor(annotation_col$NID2)
b2p1 <- colorRampPalette(c("green", "red"))
ann_colors = list(
  NID2 = c("Low" = "#2b83ba", "High" = "#d7191c")
)
pdf(file="Figure_6D.pdf",width=12,height=2)
pheatmap(leukemia_es, annotation_col = annotation_col, show_colnames = F, cluster_rows = F, cluster_cols = T, annotation_colors = ann_colors, col = b2p1(100), clustering_distance_cols = "euclidean", 
         cutree_cols = 2,
         clustering_method = "ward.D"
)
dev.off()
rownames(leukemia_es)

####
colnames(leukemia_es)
rownames(geneExp)

## start pearson correlation coeff
temp <- t(geneExp)
head(temp)
cordat<-rbind(leukemia_es,temp)
dim(cordat)

cordat[1:7, 1:3]
cordat[1:6, 1:3]
head(cordat)
cordat <- cordat[-7, ]
rownames(cordat)
batch_cor <- function(gene){
  y <- as.numeric(cordat[gene,])
  rownames <- rownames(cordat)
  do.call(rbind,future_lapply(rownames, function(x){
    dd  <- cor.test(as.numeric(cordat[x,]),y,type="spearman")
    data.frame(gene=gene,mRNAs=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}
library(future.apply)
system.time(dd <- batch_cor("NID2"))
head(dd)
nrow(dd)
dd <- dd[-nrow(dd), ]
dd$`abs(cor)`<-abs(dd$cor)
dd$`Correlation Coefficient` <- dd$cor
dd$p.value <- ifelse(dd$p.value < 0.001, "***", dd$p.value)
dd = dd[order(dd[,3]),]

dd[-ncol(dd), ]
library(ggpubr)
p<-ggdotchart(dd[-ncol(dd), ], x = "mRNAs", y = "Correlation Coefficient",
              # color = "p.value",
              label = dd[-ncol(dd), ]$p.value,
              sorting = "descending",                      
              add = "segments",                            
              add.params = list(color = "lightgray", size = 1.5),
              rotate = TRUE,                                
              dot.size = "abs(cor)",
              ggtheme = theme_pubr(),                        
              xlab=""
              
)
pdf(file="Supplementary_Figure_8C.pdf",width=12,height=3)
p <-p + gradient_color(c("blue","red")) + grids(linetype = "dashed")
p + grids(linetype = "dashed")
dev.off()

