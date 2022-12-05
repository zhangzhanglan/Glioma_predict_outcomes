library("ggpubr")
library(ggplot2)

## ## survival_cgga_Figure_3EFGHIJ.R
data <- df
data <- subset(data, Histology != "NA")
head(data)
table(data$Histology)
data$Histology <- gsub("rAOA", "Anaplastic Oligoastrocytoma", data$Histology)
data$Histology <- gsub("AOA", "Anaplastic Oligoastrocytoma", data$Histology)
data$Histology <- gsub("rOA", "Oligoastrocytoma", data$Histology)
data$Histology <- gsub("OA", "Oligoastrocytoma", data$Histology)
data$Histology <- gsub("rAA", "Anaplastic astrocytoma", data$Histology)
data$Histology <- gsub("AA", "Anaplastic astrocytoma", data$Histology)
data$Histology <- gsub("rAO", "Anaplastic oligodendroglioma", data$Histology)
data$Histology <- gsub("AO", "Anaplastic oligodendroglioma", data$Histology)
data$Histology <- gsub("rA", "Astrocytoma", data$Histology)
data$Histology <- gsub("rO", "Oligodendroglioma", data$Histology)
data$Histology <- gsub("^O$", "Oligodendroglioma", data$Histology)
data$Histology <- gsub("^A$", "Astrocytoma", data$Histology)
data$Histology <- gsub("rGBM", "Glioblastoma", data$Histology)
data$Histology <- gsub("sGBM", "Glioblastoma", data$Histology)
data$Histology <- gsub("GBM", "Glioblastoma", data$Histology)
head(data)
table(group2$IDH_mutation_status)
table(data$`1p19q_codeletion_status`)
table(group2$Label)
group2 <- data
group2 <- cbind(group2, 0)
colnames(group2)[ncol(group2)] <- "Label"
group2 <- group2[!is.na(group2$IDH_mutation_status),]
group2 <- group2[!is.na(group2$`1p19q_codeletion_status`),]
group2[((group2$Histology == "Glioblastoma") & (group2$IDH_mutation_status == "Wildtype")), ncol(group2)] = "GBM IDH-wildtype"
group2[((group2$Histology != "Glioblastoma") & (group2$IDH_mutation_status == "Mutant") & (group2$`1p19q_codeletion_status` == "Codel")), ncol(group2)] = "LGG IDH-mutant & codeleted"
group2[((group2$Histology != "Glioblastoma") & (group2$IDH_mutation_status == "Mutant") & (group2$`1p19q_codeletion_status` == "Non-codel")), ncol(group2)] = "LGG IDH-mutant & non-codeleted"
group2[((group2$Histology != "Glioblastoma") & (group2$IDH_mutation_status == "Wildtype")), ncol(group2)] = "LGG IDH-wildtype"
group2 <- subset(group2, Label != 0)
## data prepare end skip to start

## start
library(vioplot) 
head(data)
data <- group2
colnames(data)
gene = "NID2"
type = "Histology"  ## C
# type = "Label"    ## D
# type = "Grade"    ## B
colnames(data)[which(colnames(data) == gene)] <- "expression"
colnames(data)[which(colnames(data) == type)] <- "type"
table(data$type)
# seq <- c("LGG IDH-mutant & codeleted", "LGG IDH-mutant & non-codeleted", "LGG IDH-wildtype", "GBM IDH-wildtype") ## D
# seq <- c("WHO II", "WHO III", "WHO IV")    ## B
seq <- c("Oligoastrocytoma", "Oligodendroglioma", "Astrocytoma", "Glioblastoma")  ## C
data <- data[which((data$type) %in% seq), ]
ordered.item <- factor(1:length(seq),labels = seq)
data$type <- factor(data$type, levels = levels(ordered.item))
list = table(data$type)
table(data$type)
num <- as.vector(list)
num
seq_certain = seq[seq %in% data$type]
seq_certain
# n = length(unique(data$Histology))
n = length(unique(seq_certain))
dist = max(range(data$expression))/10
ylim = range(data$expression)
ylim[2] = ylim[2]+dist
ylim

##
pdf("Figure_3B.pdf",width = 5,height = 5)
pdf("Figure_3C.pdf",width = 5,height = 5)
pdf("Figure_3D.pdf",width = 5,height = 5)
# pdf(paste(gene, "_", type, ".pdf", sep =""), width = 5, height = 5)
par(mar=c(3.6, 4, 1.6, 2.1), mgp = c(3, 1, 0))
boxplot(-1, -1 , xlim=c(0.5, n + 0.5), ylim=ylim, xaxt="n", col="yellow")
seq_certain[1]
data[data$type == seq_certain[2], 3]
color = "#2b83ba"
labs <- paste(seq_certain, "\n(N=", num, ")", sep = "")
# labs <- c("LGG   \nIDH-mutant\ncodeleted\n(N=166)", "LGG   \nIDH-mutant\nnon-codeleted\n(N=246)", "LGG\nIDH-wildtype\n(N=131)", "GBM\nIDH-wildtype\n(N=261)")     ## D
for(i in 1:n){
  vioplot(x = data[data$type == seq_certain[i], "expression"], at = i, add = T, col = color)   ## NID2 col 3
}
colnames(data)
aov_model = aov(expression ~ type, data = data)
stat_result = summary(aov_model)[[1]][1,][4:5]
text(x = n + 0.5 ,y=max(data$expression)*1.05, labels = paste("F value = ",signif(stat_result[1,1],3),"\nPr(>F) = ",
    signif(stat_result[1,2],3),sep=""),cex = 1.1,col = "black",pos = 2)
axis(1, at = 1:n, labels = FALSE, tick = TRUE)
text(x = 1:length(seq_certain),
     y = par("usr")[3] - 0.7,
     labels =labs,
     xpd = NA,
     srt = 20,   ## C
     cex = 0.8)
title(ylab = paste(gene, "Expression", sep = " "))
blackhole = dev.off()
paste(gene, "Expression", sep = " ")
  
