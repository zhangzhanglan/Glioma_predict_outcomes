require(GEOquery)
require(Biobase)
#### GSE7696
eset <- getGEO("GSE7696", destdir = './',AnnotGPL = T,getGPL = F)
beta.m <- exprs(eset[[1]])
pD.all <- pData(eset[[1]])
colnames(pD.all)
pD <- pD.all[, c("title", "geo_accession", "characteristics_ch1.1")]
head(beta.m)
table(pD$characteristics_ch1.1)
a <- beta.m["204114_at", ]
a <- as.data.frame(a)
head(a)
colnames(a) <- "expression"
colnames(a) <- "NID2"
head(pD)
a
rt <- merge(a, pD, by= "row.names")
table(rt$Histology)
rt$Histology <- rt$characteristics_ch1.1
rt$Histology <- gsub("disease status: ", "", rt$Histolog)
rt$Histology <- gsub("re-recurrent GBM", "GBM", rt$Histolog)
rt$Histology <- gsub("recurrent GBM", "GBM", rt$Histolog)
rt$Histology <- gsub("non-tumoral", "non-tumor", rt$Histolog)

#### GSE4290
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)
eset <- getGEO("GSE4290", destdir = './',AnnotGPL = T,getGPL = F)
beta.m <- exprs(eset[[1]])
pD.all <- pData(eset[[1]])
colnames(pD.all)
pD <- pD.all[, c("title", "geo_accession", "Histopathological diagnostic:ch1", "characteristics_ch1")]
table(pD$`Histopathological diagnostic:ch1`)
table(pD$characteristics_ch1)
a <- beta.m["204114_at", ]
a <- as.data.frame(a)
head(a)
colnames(a) <- "expression"
colnames(a) <- "NID2"
head(pD)
a
rt <- merge(a, pD, by= "row.names")
table(rt$Histology)
rt$Histology <- rt$`Histopathological diagnostic:ch1`
rt$Histology <- gsub("astrocytoma, grade 2", "LGG", rt$Histolog)
rt$Histology <- gsub("astrocytoma, grade 3", "LGG", rt$Histolog)
rt$Histology <- gsub("oligodendroglioma, grade 2", "LGG", rt$Histolog)
rt$Histology <- gsub("oligodendroglioma, grade 3", "LGG", rt$Histolog)
rt$Histology <- gsub("glioblastoma, grade 4", "GBM", rt$Histolog)
rt <- subset(rt, Histology != "")
rt <- subset(rt, Histology != "LGG")
head(rt)
rownames(rt) <- rt$Row.names
rt[, c(2, 7)]
rt <- rt[, c(2, 7)]
library(edgeR)
rt$NID2 <- as.numeric(log2(cpm(rt$NID2)+1))


####GSE16011
require(GEOquery)
require(Biobase)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)
eset <- getGEO("GSE16011", destdir = './',AnnotGPL = T,getGPL = F)
beta.m <- exprs(eset[[1]])
boxplot(beta.m) 
pD.all <- pData(eset[[1]])
colnames(pD.all)
pD <- pD.all[, c("title", "geo_accession", "characteristics_ch1.2")]
p = identical(rownames(pD),colnames(beta.m))
if(!p) beta.m <- beta.m[,match(rownames(pD),colnames(beta.m))]
head(beta.m)
table(pD$characteristics_ch1.2)
a <- beta.m["22795_at", ]  ## grep -i nitrogen GPL8542.txt
a <- as.data.frame(a)
head(a)
colnames(a) <- "expression"
head(pD)
a
rt <- merge(a, pD, by= "row.names")
colnames(rt)
table(rt$Histology)
rt
rt$Histology <- rt$characteristics_ch1.2
rt$Histology <- gsub("^$", "Normal adult brain", rt$Histolog)
library(ggpubr)

temp_ggplot <- as.data.frame(rt)
colnames(temp_ggplot)
colnames(temp_ggplot) <- c("Row.names", "NID2", "group", "geo_accession", "characteristics_ch1.2", "Histology")
temp_ggplot$NID2 <- as.numeric(temp_ggplot$NID2 )
temp_ggplot$group <- gsub("glioma [0123456789]+", "glioma", temp_ggplot$group)
temp_ggplot$group <- gsub("control [0123456789]+", "control", temp_ggplot$group)
temp_ggplot$Histology <- gsub("histology: GBM \\(grade IV\\)", "GBM", temp_ggplot$Histology)
temp_ggplot$Histology <- gsub("histology: A \\(grade [I]+\\)", "LGG", temp_ggplot$Histology)
temp_ggplot$Histology <- gsub("histology: OA \\(grade [I]+\\)", "LGG", temp_ggplot$Histology)
temp_ggplot$Histology <- gsub("histology: PA \\(grade I\\)", "LGG", temp_ggplot$Histology)
temp_ggplot$Histology <- gsub("histology: OD \\(grade [I]+\\)", "LGG", temp_ggplot$Histology)

table(temp_ggplot$group)
compare_means(NID2 ~ group,  data = temp_ggplot, method = "t.test")
table(temp_ggplot$group)
group_list = factor(temp_ggplot$group,
                    levels = c("control","glioma"))

library(limma)
design=model.matrix(~group_list)
dim(beta.m)
design=model.matrix(~temp_ggplot$group)
fit=lmFit(beta.m, design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
library(dplyr)
deg <- mutate(deg, probe_id=rownames(deg))
deg["22795_at", ]
file <- "/Users/lanzhzh/Downloads/GPL8542.txt"
gpl <- read.table(file,sep="\t", header=T, check.names=F,row.names=1, skip = "#")
colnames(gpl)
head(deg)
deg <- merge(deg, gpl, by= "row.names")
deg <- deg[!duplicated(deg$Description),]
logFC_t = 2
P.Value_t = 0.001
k1 = (deg$P.Value < P.Value_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value < P.Value_t)&(deg$logFC > logFC_t)
change = ifelse(k1, "Down",ifelse(k2, "Up", "Stable"))
deg <- mutate(deg,change)
head(deg)
library(dplyr)
library(ggplot2)
dat  = deg
for_label <- dat%>% 
  filter(Description %in% "nidogen 2 (osteonidogen)") 
for_label$Description <- "NID2"
p <- ggplot(data = dat,
            aes(x = logFC,
                y = -log10(P.Value))) +
  geom_point(alpha=1, size=3,
             aes(color=change)) +
  xlab("log2 fold change") + ylab("-log10 P-value") +
  ylim(0, 40) + 
  scale_color_manual(values=c("#2b83ba", "black","#d7191c"))+
  theme(plot.title = element_text(size=15, hjust = 0.5))+
  theme_bw()
volcano_plot <- p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = Description), size = 5, 
    data = for_label,
    color="black"
  )
pdf(file="Figure_1B.pdf", width=5, height=5)
volcano_plot
dev.off()

#### WGCNA
dat=avereps(beta.m)
dat=dat[apply(dat,1,sd)>0.01,]
dat=as.data.frame(t(dat))  # 转换为行是样本，列是基因
library(FactoMineR)
library(factoextra) 
library("WGCNA")

datExpr0 <- dat
gsg = goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

dim(temp_ggplot)
dim(datExpr0)
samplesame <- intersect(temp_ggplot$Row.names, row.names(datExpr0))
samplesame
temp_ggplot <- temp_ggplot[temp_ggplot$Row.names %in% samplesame, ]
dat.pca <- datExpr0[samplesame, ]

## 7k pca
dat.pca <- t(sigExp)    #### sigExp from WGCNA_Figure_S1.R
dat.pca <- PCA(dat.pca, graph = FALSE,  scale = TRUE)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom = 'point',
                         col.ind = temp_ggplot$group, 
                         palette = c("#2b83ba", "#d7191c"), 
                         addEllipses = TRUE,
                         legend.title = "Groups"
)
table(temp_ggplot$group)
pdf(file="Figure_1A.pdf", width=6, height=5)
pca_plot
dev.off()
