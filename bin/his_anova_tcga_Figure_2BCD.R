library("ggpubr")
library(ggplot2)

## survival_tcga_Figure_2EFGH.R
data <- df
colnames(df)
table(data$neoplasm_histologic_grade)
table(data$primary_diagnosis.diagnoses)
data <- subset(data, primary_diagnosis.diagnoses != "")

head(data)
colnames(group2)
table(group2$ldh1_mutation_found)
table(group2$primary_diagnosis.diagnoses)
table(group2$Label)
group2 <- data
group2 <- cbind(group2, 0)
colnames(group2)[ncol(group2)] <- "Label"
# group2 <- group2[!is.na(group2$ldh1_mutation_found),]
group2 <- subset(group2, ldh1_mutation_found != "")
# group2 <- subset(group2, ldh1_mutation_found != "NA")

##
group2[((group2$primary_diagnosis.diagnoses == "Glioblastoma") & (group2$ldh1_mutation_found == "NA")), ncol(group2)] = "GBM, NOS"
group2[((group2$primary_diagnosis.diagnoses != "Glioblastoma") & (group2$ldh1_mutation_found == "YES")), ncol(group2)] = "LGG IDH-mutant"
group2[((group2$primary_diagnosis.diagnoses != "Glioblastoma") & (group2$ldh1_mutation_found == "NO")), ncol(group2)] = "LGG IDH-wildtype"
group2 <- subset(group2, Label != 0)

##
head(data)
table(data$primary_diagnosis.diagnoses)
data$primary_diagnosis.diagnoses <- gsub("Astrocytoma, NOS", "Astrocytoma", data$primary_diagnosis.diagnoses)
data$primary_diagnosis.diagnoses <- gsub("Oligodendroglioma, NOS", "Oligodendroglioma", data$primary_diagnosis.diagnoses)
data$primary_diagnosis.diagnoses <- gsub("Mixed glioma", "Oligoastrocytoma", data$primary_diagnosis.diagnoses)
table(data$primary_diagnosis.diagnoses)
## anova
library(vioplot) 
head(data)
colnames(data)
# data <- group2
colnames(data)
gene = "NID2"
type = "neoplasm_histologic_grade"    ## B
type = "primary_diagnosis.diagnoses"    ## C
type = "Label"    ## D
colnames(data)[which(colnames(data) == gene)] <- "expression"
colnames(data)[which(colnames(data) == type)] <- "type"
table(data$type)
# mean(data[data$type == "G4", ]$expression)
# mean(data[data$type == "G3", ]$expression)
# mean(data[data$type == "G2", ]$expression)
seq <- c("G2", "G3", "G4")  ## B
seq <- c("LGG IDH-mutant", "LGG IDH-wildtype", "GBM, NOS")  ## D
seq <- c("Astrocytoma", "Oligodendroglioma", "Oligoastrocytoma", "Glioblastoma")    ## C
data <- data[which((data$type) %in% seq), ]
ordered.item <- factor(1:length(seq),labels = seq)
data$type <- factor(data$type, levels = levels(ordered.item))
list = table(data$type)
table(data$type)
num <- as.vector(list)
num
seq_certain = seq[seq %in% data$type]
seq_certain
n = length(unique(seq_certain))
dist = max(range(data$expression))/10
ylim = range(data$expression)
ylim[2] = ylim[2]+dist
##
# pdf(paste(gene, "_", type, ".pdf", sep =""), width = 5, height = 5)
# pdf("Figure_2B.pdf", width = 5, height = 5)
pdf("Figure_2C.pdf", width = 5, height = 5)
pdf("Figure_2D.pdf", width = 5, height = 5)
par(mar=c(3.6, 4, 1.6, 2.1), mgp = c(3, 1, 0))
boxplot(-1,-1 , xlim=c(0.5,n + 0.5), ylim=ylim, xaxt="n", col="yellow")
seq_certain[1]
data[data$type == seq_certain[2], 3]
color = "#2b83ba"
labs <- paste(seq_certain, "\n(N=", num, ")", sep = "")
for(i in 1:n){
  vioplot(x = data[data$type == seq_certain[i], 3], at = i, add = T, col = color)   ## NID2 in col 3
}
aov_model = aov(expression ~ type, data = data)
stat_result = summary(aov_model)[[1]][1,][4:5]
text(x = n + 0.5 ,y=max(data$expression)*1.05, labels = paste("F value = ",signif(stat_result[1,1],3),"\nPr(>F) = ",
    signif(stat_result[1,2],3),sep=""),cex = 1.1,col = "black",pos = 2)
axis(1, at = 1:n, labels = FALSE, tick = TRUE)
text(x = 1:length(seq_certain),
     y = par("usr")[3] - 0.7,
     labels =labs,
     xpd = NA,
     # srt = 20,   ## B xx
     cex = 0.8)
title(ylab = paste(gene, "Expression", sep = " "))
blackhole = dev.off()
paste(gene, "Expression", sep = " ")
  
