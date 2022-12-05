exp <- read.csv("../data/TMA.csv")
colnames(exp)
dim(exp)
head(exp)
exp <- exp[exp$组织类型 != "正常大脑", ]
table(exp$intensity)
table(exp$IRS)
exp$年龄
exp$年龄 <- gsub("成年", 18, exp$年龄)
mean(as.numeric(exp$年龄))
age <- ifelse(exp$年龄 <= 40,"<=40",">40")
table(age)
table(exp$病理分型)

table(exp$病理分级)
table(exp$病理分型)
exp[exp$病理分级=="Ⅲ-Ⅳ级", ]
exp[exp$病理分级=="Ⅱ-Ⅲ级", ]
exp[exp$蜡块编号 == "P01A0445-B30-C01", "病理分级"] <- "Ⅱ级"
exp[exp$蜡块编号 == "P01A0444-B30-C01", "病理分级"] <- "Ⅱ级"
exp[exp$病理分型 == "间变型星形细胞瘤+胶质母细胞瘤", ]
exp$who <- ifelse(exp$组织类型 == "正常大脑", "Normal", exp$病理分级)
table(exp$type)
exp$type <- exp$who
exp$type <- gsub("Ⅱ-Ⅲ级", "LGG", exp$type)
exp$type <- gsub("Ⅲ-Ⅳ级", "LGG", exp$type)  
exp$type <- gsub("Ⅲ级", "LGG", exp$type)
exp$type <- gsub("Ⅱ级", "LGG", exp$type)
exp$type <- gsub("Ⅰ级", "LGG", exp$type)
exp$type <- gsub("Ⅳ级", "GBM", exp$type)
exp$who <- gsub("Ⅱ-Ⅲ级", "WHO III", exp$who)
exp$who <- gsub("Ⅲ-Ⅳ级", "WHO III", exp$who)
exp$who <- gsub("Ⅲ级", "WHO III", exp$who)
exp$who <- gsub("Ⅱ级", "WHO II", exp$who)
exp$who <- gsub("Ⅰ级", "WHO I", exp$who)
exp$who <- gsub("Ⅳ级", "WHO IV", exp$who)
exp$组织类型
exp$sex <- exp$性别
exp$sex <- ifelse(exp$性别 == "男", "M", "F")
table(exp$sex)
mean(as.numeric(exp$年龄))
round(sd(as.numeric(exp$年龄)), 2)
list <- c("GBM", "LGG", "Normal")
outTab=data.frame()
for (i in list) {
  mean <- round(mean(as.numeric(exp[exp$type== i, ]$年龄)), 2)
  sd <- round(sd(as.numeric(exp[exp$type== i, ]$年龄)), 2)
  outTab <- rbind(outTab,cbind(Type=i, mean, sd))
}
outTab

wilcox.test(as.numeric(exp[exp$type== "GBM", ]$年龄), as.numeric(exp[exp$type== "LGG", ]$年龄))$p.value
wilcox.test(as.numeric(exp[exp$type== "GBM", ]$年龄), as.numeric(exp[exp$type== "LGG", ]$年龄))
as.numeric(exp[exp$type== "GBM", ]$年龄)
as.numeric(exp[exp$type== "LGG", ]$年龄)     
exp$年龄
exp$who
exp$NID2
exp$HBraG125PG02数据PDL1
exp$IRS
exp$Final.IRS.score..0.12.
colnames(exp)
# data <- data.frame(as.numeric(exp$年龄), exp$type, exp$who, exp$sex, exp$IRS, exp$NID2, exp$PDL1Final, exp$PDL.IRS, exp$病理分型)   ## heatmap
data <- data.frame(as.numeric(exp$年龄), exp$type, exp$who, exp$sex, exp$Final.IRS.score..0.12., exp$NID2, exp$PDL1Final, exp$PDL.IRS, exp$病理分型)   ## IRS plot
colnames(data) <- c("age", "type", "who", "sex", "IRS",  "NID2", "PDL1_s", "PDL1", "histological")
data$NID2_g <- data$NID2
data$NID2_g <- ifelse(data$NID2 < 2, "Low NID2 expression", "High NID2 expression")
table(data$PDL1)
exp[exp$IRS == "mild", ]
exp[exp$蜡块编号 == "P01A0277-B30-C01", ]
table(exp$IRS)
data <- subset(data, type != "Normal")
table(data$histological)
data$histological <- gsub("GBM", "Glioblastoma", data$histological)
# data$histological <- gsub("间变型星形细胞瘤+胶质母细胞瘤", "Anaplastic astrocytoma", data$histological)
data$histological <- gsub("节细胞胶质瘤", "Gangliogliomas", data$histological)
data$histological <- gsub("胶质母细胞瘤，胶质肉瘤", "Glioblastoma", data$histological)
data$histological <- gsub("胶质母细胞瘤", "Glioblastoma", data$histological)
data$histological <- gsub("间变型少突胶质瘤", "Anaplastic oligodendroglioma", data$histological)
data$histological <- gsub("间变型星形细胞瘤", "Anaplastic astrocytoma", data$histological)
data$histological <- gsub("少突胶质瘤", "Oligodendroglioma", data$histological)
data$histological <- gsub("间变型少突星形细胞瘤", "Oligo astrocytoma", data$histological)
data$histological <- gsub("少突星形细胞瘤，部分间变型", "Oligo astrocytoma", data$histological)
data$histological <- gsub("少突星形细胞瘤", "Oligo astrocytoma", data$histological)
data$histological <- gsub("毛细胞型星形细胞瘤", "Pilocytic astrocytoma", data$histological)
data$histological <- gsub("多形性黄色星形细胞瘤", "Pleomorphic xanthoastrocytoma", data$histological)
data$histological <- gsub("弥漫型星形细胞瘤", "Diffuse astrocytoma", data$histological)
table(data$sex)
table(data$IRS)
table(data$histological)
data[data$histological == "Anaplastic astrocytoma+Glioblastoma", ]$histological <- "Mixed glioma"
data$histological <- gsub("Diffuse astrocytoma", "Astrocytoma, NOS", data$histological)
data$histological <- gsub("Pleomorphic xanthoastrocytoma", "Astrocytoma, NOS", data$histological)
data$histological <- gsub("Pilocytic astrocytoma", "Astrocytoma, NOS", data$histological)
data$histological <- gsub("Anaplastic astrocytoma", "Anaplastic, astrocytoma", data$histological)
data$histological <- gsub("Oligodendroglioma", "Oligodendroglioma, NOS", data$histological)
data$histological <- gsub("Anaplastic oligodendroglioma", "Oligodendroglioma, anaplastic", data$histological)
data$histological <- gsub("Gangliogliomas", "Mixed glioma", data$histological)
data$histological <- gsub("Oligo astrocytoma", "Mixed glioma", data$histological)
table(data$histological)

### data ready
# write.table(data,file="TMA_u.txt",sep="\t",quote=F,col.names=F)       

wilcox.test(age ~ type, data)
par(mfrow=c(1,2))
hist(data[data$type =='LGG', ]$age, ylab = 'Frequency for Age', main = "",
     xlab=paste('LGG ', paste(round(mean(as.numeric(exp[exp$type== "LGG", ]$年龄)), 2),
                round(sd(as.numeric(exp[exp$type== "LGG", ]$年龄)), 2), sep = "±")))
hist(data[data$type =='GBM', ]$age, ylab = 'Frequency for Age', main = "",
     xlab=paste('GBM ', paste(round(mean(as.numeric(exp[exp$type== "GBM", ]$年龄)), 2),
                              round(sd(as.numeric(exp[exp$type== "GBM", ]$年龄)), 2), sep = "±")))
library(ggplot2)
ggplot(data,aes(x=type)) + geom_bar(aes(fill=factor(sex)),position="dodge") + xlab("") 
data$type
## 7 x 5 in

table=matrix(c(sum(as.numeric(exp[exp$type== "LGG", ]$年龄) > mean(as.numeric(exp$年龄))),
               sum(as.numeric(exp[exp$type== "GBM", ]$年龄) > mean(as.numeric(exp$年龄))), 
               sum(as.numeric(exp[exp$type== "LGG", ]$年龄) < mean(as.numeric(exp$年龄))), 
               sum(as.numeric(exp[exp$type== "GBM", ]$年龄) < mean(as.numeric(exp$年龄)))),nrow=2)
colnames(table)=c("High","Low")
rownames(table)=c("LGG","GBM")
table
fisher.test(table)$p.value

table <- xtabs(~sex+type, data=data)
table <- xtabs(~NID2+type, data=data)
table <- xtabs(~NID2_g + type, data = data)
table
table=matrix(c(10, 9, 17, 33), nrow = 2)
table=matrix(c(7, 11, 11, 33), nrow = 2)
table=matrix(c(8, 5, 10, 33), nrow = 2)
fisher.test(table, alternative = "two.sided")$p.value
table
ftable(table)
summary(table)
prop.table(table, 2)

table(data$type)
colnames(data)
as.numeric(data[data$type== "GBM", ]$NID2)
pval <- wilcox.test(as.numeric(data[data$type== "GBM", ]$NID2), as.numeric(data[data$type== "LGG", ]$NID2))$p.value
if(pval<0.001){
  pval="< 0.001"
}else{
  pval=paste0("=",sprintf("%.03f",pval))
}
colnames(data)
data[1:5, 1:9]
as.numeric(data[data$type== "GBM", ]$age)
as.numeric(data[data$type== "LGG", ]$IRS)
colnames(data)
GBM_N_s <- t.test(as.numeric(data[data$type== "GBM", ]$IRS))
LGG_N_s <- t.test(as.numeric(data[data$type== "LGG", ]$IRS))
df <- matrix(c("GBM", "LGG", GBM_N_s$estimate, LGG_N_s$estimate, GBM_N_s$conf.int[1], LGG_N_s$conf.int[1],  GBM_N_s$conf.int[2], LGG_N_s$conf.int[2]), nrow = 2)
df
colnames(df) <- c("type", "mean", "sd1", "sd2")
df <- as.data.frame(df)
df$type <- as.factor(df$type)
df$mean <- as.numeric(df$mean)
df$mean
df$sd1 <- as.numeric(df$sd1)
df$sd2
df$sd2 <- as.numeric(df$sd2)

library(ggplot2)   ## final
pdf(file="Figure_4E.pdf", width=5, height=5)
table(data$type)
df$type <- c("GBM (N=53)", "LGG (N=67)")
seq <- c("LGG (N=67)", "GBM (N=53)")
ordered.item <- factor(1:length(seq),labels = seq)
df$type <- factor(df$type, levels = levels(ordered.item))
table(data$type)
df
p<- ggplot(df, aes(x=type, y=mean, fill=type)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=sd1, ymax=sd2), width=.2, position=position_dodge(.9)) 
p+ylim(0, 8) + labs(title="", x="", y = "NID2 Immunoreactive Score")+
  theme_classic() +
  theme(legend.position="none") +
  scale_fill_manual(values=c('#2b83ba','#d7191c')) +
  annotate("text", label = paste("Wilcox.test, P-value ", pval, sep=""), x= 1.5 , y=8, cex = 5)
dev.off()
