library(ComplexHeatmap)
## bar_TMA_Figure_4E.R
head(data)
df <- data
df <- subset(df, histological != "NA")
df
df$who <- gsub("WHO ", "", df$who)
df <- df[order(df$who, df$NID2),]
colnames(df)
Risk <- data.frame(data$age, data$type, data$who, data$sex, data$histological, data$PDL1, data$NID2)

table(data$NID2)
df$IRS
head(Risk)
row.names(Risk)
colnames(Risk) <- c("Age", "Subtype", "Grade", "Gender", "Histological", "PD-L1", "NID2")
table(Risk$Grade)
Risk$Grade <- gsub("WHO ", "", Risk$Grade)
table(Risk$Grade)
colnames(Risk)
Risk <- Risk[order(Risk$Grade, Risk$NID2),]
rt <- Risk
sigVec=c()
head(rt)
rt
for(clinical in colnames(rt[,1:(ncol(rt)-1)])){
	data=rt[c("NID2", clinical)]
	colnames(data)=c("NID2", "clinical")
	data=data[(data[,"clinical"]!="unknow"),]
	tableStat=table(data)
	stat=chisq.test(tableStat, simulate.p.value = TRUE)
	pvalue=stat$p.value
	Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
	sigVec=c(sigVec, paste0(clinical, Sig))
}
sigVec=c(sigVec, "NID2")
sigVec
colnames(rt)=sigVec
rt$NID2 <- df$IRS
rt$NID2=factor(rt$NID2, levels=c("negative", "mild", "moderate", "strongly postive"))
rt$`PD-L1***`=factor(rt$`PD-L1***`, levels=c("negative","mild", "moderate", "strongly postive"))
rt$`Subtype***` = factor(rt$`Subtype***`, levels = c("LGG", "GBM"))
table(rt$`Histological***`)
rt$`Histological***` = factor(rt$`Histological***`, levels = c("Astrocytoma, NOS", "Anaplastic, astrocytoma", "Oligodendroglioma, NOS", "Oligodendroglioma, anaplastic", "Mixed glioma", "Glioblastoma"))

bioCol=c("#0066FF","#FF9900","#FF0000","#ed1299", "#0dbc21", "#246b93", "#cc8e12", "#d561dd", "#c93f00", 
         "#ce2523", "#f7aa5d", "#9ed84e", "#39ba30", "#6ad157", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
         "#1a918f", "#7149af", "#ff66fc", "#2927c4", "#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
         "#4aef7b", "#e86502",  "#99db27", "#e07233", "#8249aa","#cebb10", "#03827f", "#931635", "#ff523f",
         "#edd05e", "#6f25e8", "#0dbc21", "#167275", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
         "#dd27ce", "#07a301", "#ddd53e",  "#391c82", "#2baeb5","#925bea", "#09f9f5",  "#63ff4f")
colorList=list()
j=0
sigVec
library(RColorBrewer)
brewer.pal(11, "PiYG")
colorList[["Gender"]] <- c("F" = "#f1b6da", "M" = "#b8e186")
colorList[["Subtype***"]] <- c("LGG" = "#2b83ba", "GBM" = "#d7191c")
colorList[["NID2"]]=c("negative"="#abdda4", "mild"="#fdae61", "moderate" = "#2b83ba", "strongly postive" = "#d7191c")
colorList[["PD-L1***"]]=c("negative"="#abdda4", "mild"="#fdae61", "moderate" = "#2b83ba", "strongly postive" = "#d7191c")
colorList[["Histological***"]]=c("Astrocytoma, NOS"="#1a9850", "Anaplastic, astrocytoma"="#66bd63", "Oligodendroglioma, NOS"="#a6d96a", "Oligodendroglioma, anaplastic"="#d9ef8b", "Mixed glioma"="#fee08b", "Glioblastoma"="#d73027")
colorList[["Grade***"]]=c("I" = "#1a9850", "II" = "#91cf60", "II-III" = "#d9ef8b", "III" = "#fee08b", "III-IV" = "#fc8d59", "IV" = "#d73027")
table(rt$`Histological***`)
ha=HeatmapAnnotation(df=rt, col=colorList)
zero_row_mat=matrix(nrow=0, ncol=nrow(rt))
Hm=Heatmap(zero_row_mat, top_annotation=ha)
rt

pdf(file="Figure_4F.pdf", width=10, height=5)
draw(Hm, merge_legend=TRUE, heatmap_legend_side="bottom", annotation_legend_side="bottom")
dev.off()

