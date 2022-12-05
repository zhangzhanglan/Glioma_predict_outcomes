## http://www.cgga.org.cn/download.jsp
clinic1<-read.table("CGGA.mRNAseq_325_clinical.20200506.txt",header = T,check.names = F,sep = "\t")
clinic3<-read.table("CGGA.Methyl_array_159_clinical.20200506.txt",header = T,check.names = F,sep = "\t")
# clinic4<-read.table("CGGA.WEseq_286_clinical.20200506.txt",header = T,check.names = F,sep = "\t")
clinic2<-read.table("CGGA.mRNAseq_693_clinical.20200506.txt",header = T,check.names = F,sep = "\t")
clinic <- rbind(clinic1, clinic2)
samname <- intersect(clinic$CGGA_ID, clinic3$CGGA_ID)

# get exp from survival_cgga_Figure_3EFGHIJ.R
exp[1:2, 1:2]
ins
file3 = "CGGA.Methyl_array_159_gene_level.20200506.txt"
rt3=read.table(file3,sep="\t",header=T,check.names=F)
rt3=as.matrix(rt3)
rownames(rt3)=rt3[,1]
exp3=rt3[,2:ncol(rt3)]
dimnames3=list(rownames(exp3),colnames(exp3))
data3=matrix(as.numeric(as.matrix(exp3)),nrow=nrow(exp3),dimnames=dimnames3)
met <-data3[row.names(data3) %in% ins, ]
met[1:2, 1:2]
rownames(met)
met
rownames(exp) <- c("NID1 Expression", "NID2 Expression")
rownames(met) <- c("NID1 Methylation", "NID2 Methylation")
exp <- t(exp)
met <- t(met)
rownames(clinic) <- clinic$CGGA_ID
clinic <- clinic[, -1]
exp = exp[samname,]
met = met[samname,]
clinic =clinic[samname,]
tem <- merge(exp, met, by= "row.names", all.x=TRUE)
rownames(tem) <- tem$Row.names
tem <- tem[, -1]
rt <- merge(tem, clinic, by= "row.names", all.x=TRUE)
rownames(rt) <- rt$Row.names
rt <- rt[, -1]
rt <- rt[, -1]
rt <- rt[, -2]
colnames(rt)
rt <- data.frame(rt$`NID2 Expression`, rt$`NID2 Methylation`, rt$Histology, rt$Grade, rt$IDH_mutation_status, rt$`1p19q_codeletion_status`, rt$MGMTp_methylation_status)
colnames(rt)
colnames(rt) <- c("NID2 Expression", "NID2 Methylation", "Group", "WHO Grade", "IDH Mutation", "1p/19q Codeletion", "MGMT Methylation")
rt$Group <- gsub("AA", "LGG", rt$Group)
rt$Group <- gsub("AO", "LGG", rt$Group)
rt$Group <- gsub("rAA", "LGG", rt$Group)
rt$Group <- gsub("rAO", "LGG", rt$Group)
rt$Group <- gsub("rA", "LGG", rt$Group)
rt$Group <- gsub("rO", "LGG", rt$Group)
rt$Group <- gsub("O", "LGG", rt$Group)
rt$Group <- gsub("A", "LGG", rt$Group)
rt$Group <- gsub("rLGG", "LGG", rt$Group)
rt$Group <- gsub("rGBM", "GBM", rt$Group)
rt$Group <- gsub("sGBM", "GBM", rt$Group)

rt <- rt[order(rt$`NID2 Expression`), ]
# rt <- rt[order(rt$`IDH Mutation`), ]
table(rt$`WHO Grade`)
# rt <- rt[order(rt$Group), ]
# rt$Group = factor(rt$Group, levels = c("LGG", "GBM"))
rt <- rt[order(rt$`WHO Grade`, decreasing = FALSE), ]
rt$`WHO Grade` <- factor(rt$`WHO Grade`)
rt$`WHO Grade`
colnames(rt) <- c("NID2 Expression", "NID2 Methylation", "Group", "WHO Grade", "IDH Mutation", "1p/19q Codeletion", "MGMT Methylation")
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")
bioCol=c("#2b83ba","#d7191c", "#91cf60", "#fee08b", "#d73027", "#fdae61","#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#c93f00", 
         "#ce2523", "#f7aa5d", "#9ed84e", "#39ba30", "#6ad157", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
         "#1a918f", "#7149af", "#ff66fc", "#2927c4", "#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
         "#4aef7b", "#e86502",  "#99db27", "#e07233", "#8249aa","#cebb10", "#03827f", "#931635", "#ff523f",
         "#edd05e", "#6f25e8", "#0dbc21", "#167275", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
         "#dd27ce", "#07a301", "#ddd53e",  "#391c82", "#2baeb5","#925bea", "#09f9f5",  "#63ff4f")
colorList=list()
j=0
for(cli in colnames(rt[,3:(ncol(rt))])){
  cliLength=length(levels(factor(rt[,cli])))
  cliCol=bioCol[(j+1):(j+cliLength)]
  j=j+cliLength
  names(cliCol)=levels(factor(rt[,cli]))
  cliCol["unknow"]="white"
  colorList[[cli]]=cliCol
}
col_fun = colorRamp2(breaks = seq(0, 1, length.out = 12),
                     colors = brewer.pal(12,"Set3"))

colorList[['NID2 Expression']] <- colorRamp2(c(min(rt$`NID2 Expression`), max(rt$`NID2 Expression`)), c("white", "red"))
colorList[['NID2 Methylation']] <- colorRamp2(c(min(rt$`NID2 Methylation`), max(rt$`NID2 Methylation`)), c("white", "red"))
rt$`NID2 Expression`
df <- rt[, -(1:3)]
df
ha=HeatmapAnnotation(df=df, col=colorList,
                     
                     "NID2 Expression" = anno_points(rt$`NID2 Expression`, height = unit(2, "cm"),
                      gp = gpar(col = "#d7191c", fill = "#d7191c")),
                     "NID2 Methylation" = anno_barplot(rt$`NID2 Methylation`, height = unit(1, "cm"),
                                                       gp = gpar(col = "#d7191c", fill = "#d7191c")),
                     border = TRUE, na_col = "black", annotation_name_side = "right", simple_anno_size_adjust = TRUE)

zero_row_mat=matrix(nrow=0, ncol=nrow(rt))
Hm=Heatmap(zero_row_mat, top_annotation=ha, width = unit(18, "cm"), column_split = c(rep(" LGG", 19), rep("GBM", 7)))
pdf(file="Supplementary_Figure_4.pdf", width=10, height=5)
table(rt$Group)
draw(Hm, merge_legend=TRUE, heatmap_legend_side="bottom", annotation_legend_side="bottom")
dev.off()
 
