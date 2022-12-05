library(DESeq2)
library(org.Hs.eg.db)
library(tibble)
library(dplyr)
library(tidyr)
library(fgsea)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)

## start begin
outTab <- read.table("../data/all_bestcutoff.txt", header = TRUE)
ranks <- as.numeric(outTab$logFC)
head(ranks)
names(ranks) <- outTab$gene
head(ranks)

## http://www.gsea-msigdb.org/gsea/login.jsp
gmtFile="c5.all.v7.5.1.symbols.gmt"
pathways.hallmark <- gmtPathways(gmtFile)
fgseaRes <- fgsea(pathways=pathways.hallmark, 
                  ranks, 
                  # nPermSimple = 10000,
                  minSize = 15, 
                  maxSize = 500)
# save(fgseaRes, file = "../data/gesa_tcga_house.Rdata")
load("../data/gesa_tcga_house.Rdata")
fgseaRes

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) 
head(outTab)
res <- outTab
colnames(res)[1] <- "SYMBOL"
head(pathways.hallmark$GOBP_ESTABLISHMENT_OF_MITOTIC_SPINDLE_ORIENTATION)
gene.in.pathway <- pathways.hallmark %>%
  enframe("pathway", "SYMBOL") %>% 
  unnest(cols = c(SYMBOL)) %>% 
  inner_join(res, by="SYMBOL", copy = TRUE)
head(gene.in.pathway)

fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
head(fgseaResTidy)
fgseaResTidy.all <- fgseaResTidy
fgseaResTidy.all$NES
head(fgseaResTidy.all)
fgseaResTidy.all[fgseaResTidy.all$pathway == "GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT", ]
fgseaResTidy.all[fgseaResTidy.all$pathway == "GOBP_REGULATION_OF_VASCULATURE_DEVELOPMENT", ]
fgseaResTidy.all[fgseaResTidy.all$pathway == "GOBP_RESPONSE_TO_BACTERIUM", ]
fgseaResTidy.all[fgseaResTidy.all$pathway == "GOBP_BLOOD_VESSEL_MORPHOGENESIS", ]
fgseaResTidy.all[fgseaResTidy.all$pathway == "GOCC_EXTERNAL_ENCAPSULATING_STRUCTURE", ]
fgseaResTidy <- fgseaResTidy.all[abs(fgseaResTidy.all$NES) > 1.8 & fgseaResTidy.all$size > 335, ]
dim(fgseaResTidy)
head(reorder(fgseaResTidy.all$pathway, fgseaResTidy.all$NES))
reorder(fgseaResTidy$pathway, fgseaResTidy$NES)
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES, fill = size)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  )
fgseaResTidy$pathway

hall
cluster <- c("M5: Transmenbrane transporter", "M3: Vessel morphogensis", "M1: Cell cycle regulation", "M2: Immune response", "M2: Immune response", "M2: Immune response", "M2: Immune response", "M1: Cell cycle regulation", "M1: Cell cycle regulation", "M2: Immune response", "M2: Immune response", "M2: Immune response", "M2: Immune response", "M5: Transmenbrane transporter", "M5: Transmenbrane transporter", "M4: External encapsulating", "M5: Transmenbrane transporter", "M1: Cell cycle regulation", "M5: Transmenbrane transporter", "M5: Transmenbrane transporter")
cluse <- cbind(pathway=hall, Module = cluster, no = seq(1, 20, by=1))
head(cluse)
fgseaResTidy
temp <- merge(fgseaResTidy, cluse, by = "pathway")
head(temp)
temp$pathway <- paste(temp$no, temp$pathway, sep = ".")
RColorBrewer::brewer.pal(8, "Set3")[2:6]
cols <- RColorBrewer::brewer.pal(8, "Set3")[2:6]
cols <- c("#3cb346", "#d75427", "#1F78B4", "#942d8d", "#eeb401")
pdf("Figure_5A_r.pdf", width=14, height=7)
# pdf("Figure_5A_r.pdf", width=35, height=7)   ## bottom legend
ggplot(temp, aes(reorder(pathway, NES), NES, fill = Module)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.text.y = element_text(size = 10), 
        legend.position= "none",
        # legend.position = "top",
        legend.title = element_text(size=30), legend.text = element_text(size=30),
        panel.grid=element_blank()) +
  coord_flip() +
  labs(x="", y="",
  ) 
dev.off()

fgseaResTidy$pathway
plotEnrichment(pathway = pathways.hallmark[["GOBP_VASCULATURE_DEVELOPMENT"]], ranks)
plotEnrichment(pathway = pathways.hallmark[["GOCC_CHROMOSOMAL_REGION"]], ranks)
plotEnrichment(pathway = pathways.hallmark[["GOCC_POSTSYNAPTIC_MEMBRANE"]], ranks)
plotEnrichment(pathway = pathways.hallmark[["GOBP_BLOOD_VESSEL_MORPHOGENESIS"]], ranks)
i = "GOBP_VASCULATURE_DEVELOPMENT"
i = "GOCC_EXTERNAL_ENCAPSULATING_STRUCTURE"
plotEnrichment(pathway = pathways.hallmark[[i]], 
               gseaParam = 1, ticksSize = 0.5, stats= ranks) + 
  labs(title=i) + theme(plot.title = element_text(hjust = 0.5, face="bold"))

plotGseaTable(pathways.hallmark[fgseaResTidy$pathway], ranks, fgseaRes, 
              gseaParam=0.5)

fgseaResTidy
library(enrichplot)

sig.path <- fgseaResTidy$pathway[fgseaResTidy$adjPvalue == "significant"]
dim(gene.in.pathway)
sig.gen <- unique(na.omit(gene.in.pathway$SYMBOL[gene.in.pathway$pathway %in% sig.path]))
sig.gen
head(gene.in.pathway)
dim(gene.in.pathway.ins)
fgseaResTidy$pathway
gene.in.pathway.ins <- gene.in.pathway[gene.in.pathway$pathway %in% fgseaResTidy$pathway, ]
gene.in.pathway.ins[gene.in.pathway.ins$pathway == "GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT", ]
### heatmap
h.dat <- dcast(gene.in.pathway.ins[, c(1,2)], SYMBOL~pathway)
h.dat[h.dat$SYMBOL == "NID2", ]
names <- h.dat$SYMBOL
rownames(h.dat)
rownames(h.dat) <- h.dat$SYMBOL
h.dat <- h.dat[, -1]
h.dat <- !is.na(h.dat)
h.dat <- (apply(h.dat, 2, as.numeric))
rownames(h.dat) <- names

h.dat <- h.dat[rownames(h.dat) %in% sig.gen, ]
h.dat <- h.dat[, colnames(h.dat) %in% sig.path]
head(h.dat)

table(data.frame(rowSums(h.dat)))
rownames(h.dat)
h.dat["NID2", ]

outDiff <- read.table("diff.xls", header = TRUE)
outDiff[outDiff$gene == "NID2", ]  
samgene <- intersect(rownames(h.dat), outDiff$gene)
samgene == "NID2"
h.dat <- h.dat[samgene, ]
head(h.dat)
dim(h.dat)
res[res$SYMBOL == "NID2", ]
topTable <- res[res$SYMBOL %in% rownames(h.dat), ]
topTable[topTable$SYMBOL == "NID2", ]
rownames(topTable) <- topTable$SYMBOL
topTableAligned <- topTable[which(rownames(topTable) %in% rownames(h.dat)),]
topTableAligned <- topTableAligned[match(rownames(h.dat), rownames(topTableAligned)),]
all(rownames(topTableAligned) == rownames(h.dat))

colnames(topTableAligned) <- c("SYMBOL", "conMean", "treatMean","log2FoldChange", "padj", "fdr" )
topTableAligned$padj <- as.numeric(topTableAligned$padj)
dfMinusLog10FDRGenes <- data.frame(-log10(
  topTableAligned[which(rownames(topTableAligned) %in% rownames(h.dat)), 'padj']))
dfMinusLog10FDRGenes[dfMinusLog10FDRGenes == 'Inf'] <- 0
dfFoldChangeGenes <- data.frame(
  topTableAligned[which(rownames(topTableAligned) %in% rownames(h.dat)), 'log2FoldChange'])
dfGeneAnno <- data.frame(dfMinusLog10FDRGenes, dfFoldChangeGenes)
colnames(dfGeneAnno) <- c('Gene score', 'Log2FC')
dfGeneAnno[,2] <- ifelse(dfGeneAnno$Log2FC > 0, 'Up-regulated',
                         ifelse(dfGeneAnno$Log2FC < 0, 'Down-regulated', 'Unchanged'))
colours <- list(
  'Log2FC' = c('Up-regulated' = 'royalblue', 'Down-regulated' = 'yellow'))
haGenes <- rowAnnotation(
  df = dfGeneAnno,
  col = colours,
  width = unit(1,'cm'),
  annotation_name_side = 'top')


dfEnrichment <- fgseaRes[, c("pathway", "NES")]
dfEnrichment <- dfEnrichment[dfEnrichment$pathway %in% colnames(h.dat)]
dd <- dfEnrichment$pathway
dd
dfEnrichment <- dfEnrichment[, -1]
rownames(dfEnrichment) <- dd
colnames(dfEnrichment) <- 'Normalized\n Enrichment score'
h.dat[data.frame(rowSums(h.dat))  == 0, ]
haTerms <- HeatmapAnnotation(
  df = dfEnrichment,
  Term = anno_text(
    colnames(h.dat),
    rot = 45,
    just = 'right',
    gp = gpar(fontsize = 8)),
  annotation_height = unit.c(unit(1, 'cm'), unit(8, 'cm')),
  annotation_name_side = 'left')
colnames(h.dat)
hmapGSEA <- Heatmap(h.dat,
                    name = 'GSEA hallmark pathways enrichment',
                    split = dfGeneAnno[,2],
                    col = c('0' = 'white', '1' = 'forestgreen'),
                    rect_gp = gpar(col = 'grey85'),
                    cluster_rows = TRUE,
                    show_row_dend = TRUE,
                    row_title = 'Top Genes',
                    row_title_side = 'left',
                    row_title_gp = gpar(fontsize = 11, fontface = 'bold'),
                    row_title_rot = 90,
                    show_row_names = TRUE,
                    row_names_gp = gpar(fontsize = 11, fontface = 'bold'),
                    row_names_side = 'left',
                    row_dend_width = unit(35, 'mm'),
                    cluster_columns = TRUE,
                    show_column_dend = TRUE,
                    column_title = 'Enriched terms',
                    column_title_side = 'top',
                    column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                    column_title_rot = 0,
                    show_column_names = FALSE,
                    show_heatmap_legend = FALSE,
                    clustering_distance_columns = 'euclidean',
                    clustering_method_columns = 'ward.D2',
                    clustering_distance_rows = 'euclidean',
                    clustering_method_rows = 'ward.D2',
                    bottom_annotation = haTerms)
pdf("GSEA_enrichment_2_tcga_house.pdf", width=13, height=22)
draw(hmapGSEA + haGenes,
     heatmap_legend_side = 'right',
     annotation_legend_side = 'right')
dev.off()


