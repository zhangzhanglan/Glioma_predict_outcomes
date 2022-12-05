library(sva)
library(limma)
# file1="../data/tcgaSymbol_LGG.txt"
# file2="../data/tcgaSymbol_GBM.txt"
# 
# rt1=read.table(file1,sep="\t",header=T,check.names=F)
# rt1=as.matrix(rt1)
# rownames(rt1)=rt1[,1]
# exp1=rt1[,2:ncol(rt1)]
# dimnames1=list(rownames(exp1),colnames(exp1))
# data1=matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames1)
# data1=avereps(data1)
load("../data/TCGA-LGG.htseq_fpkmtotpm.Rdata")   ##
data1=avereps(stad_tpm)
group=sapply(strsplit(colnames(data1),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
table(group)
data1=data1[,group==0]
colnames(data1) <- gsub("-0[12][ABC]$", "", colnames(data1))
data1[1:4, 1:4]

# rt2=read.table(file2,sep="\t",header=T,check.names=F)
# rt2=as.matrix(rt2)
# rownames(rt2)=rt2[,1]
# exp2=rt2[,2:ncol(rt2)]
# dimnames2=list(rownames(exp2),colnames(exp2))
# data2=matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)
# data2=avereps(data2)
load("../data/TCGA-GBM.htseq_fpkmtotpm.Rdata")   ##
data2=avereps(stad_tpm)
group=sapply(strsplit(colnames(data2),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
table(group)
colnames(data2[,group==1])
data2=data2[,group==0]
colnames(data2) <- gsub("-0[12][ABC]$", "", colnames(data2))
# colnames(data2) <- gsub("-11A$", "", colnames(data2))
# colnames(data2) <- gsub("[ABC]$", "", colnames(data2))
data2[1:4, 1:4]

dim(data1)
dim(data2)

sameGene=intersect(row.names(data1),row.names(data2))
data=cbind(data1[sameGene,],data2[sameGene,])

# data <- data2       ## GBM
# data <- data1       ## LGG

data <- log2(data + 1)    ## pearson survival his_anova
# data <- log10(data + 1)   ##  drug
dim(data)

table(duplicated(colnames(data)))
data <- t(data)
data <- avereps(data)
data <- t(data)

dim(data) 

ins <- c("NID1", "NID2")
exp <-data[row.names(data) %in% ins, ]
exp[1:2, 1:2]

### survival
load("../data/survival_tcga_xena.Rdata")
clinic <- phe
rownames(clinic)
colnames(clinic)
dim(clinic)
colnames(clinic)
head(clinic)
# clinic<-clinic[!is.na(clinic$OS),]
# clinic<-clinic[!is.na(clinic$OS_status),]
# clinic<-clinic[!clinic$OS<30,]
colnames(clinic) <- c("fustat", "futime")
### clinic end

load("../data/clinical_tcga_xena.Rdata")
head(phe)
age <- phe$age_at_initial_pathologic_diagnosis
age <- na.omit(age)
mean(age)
clinic <- phe
# table(phe$ldh1_mutation_found)
# phe["TCGA-QH-A65S", ]
# clinic <- clinic[rownames(phe[phe$disease_code == "LGG" & phe$ldh1_mutation_found == "YES", ]), ]
# clinic <- clinic[rownames(phe[phe$disease_code == "LGG" & phe$ldh1_mutation_found == "NO", ]), ]

rt <- exp
rt <- t(rt)
rownames(rt)
rownames(clinic)
head(rt)
head(clinic)
dim(clinic)

# library(tidyverse)
# df <- list(rt, clinic) %>% 
#   map(~ .x %>% 
#         as.data.frame %>%
#         rownames_to_column('rn')) %>% 
#   reduce(left_join, by = 'rn') %>%
#   column_to_rownames('rn')
# dim(df)
sample <- intersect(rownames(rt), rownames(clinic))
table(rownames(rt) %in% sample)

head(sample)
rt = rt[sample,]
clinic =clinic[sample,]
df <- merge(rt, clinic, by= "row.names", all.x=TRUE)
df <- as.data.frame(df)
rownames(df) <- df$Row.names
head(df)
dim(df)
df_cli <- df
df$fustat = ifelse(df$fustat == "Dead", "1", "0")   ## ROC final
df_sur <- df
### end

##### 
gene <- "NID2"
head(df)
dim(df)

# df<-df[!df$Row.names == "TCGA-06-0125",]

dim(df)
highcutoff <- 50
lowcutoff <- 50
df <- as.data.frame(df)
rownames(df) <- df$Row.names
high_df = rownames(df[df[,gene] > quantile(df[,gene],as.numeric(highcutoff)/100),1,drop = F])
df
df[df$futime == 1448, ]
df[df$futime == 1448, ]$Row.names

quantile(df[,gene], as.numeric(highcutoff)/100)
low_df = rownames(df[df[,gene] < quantile(df[,gene], as.numeric(lowcutoff)/100),1,drop = F])

high_num = length(high_df)
low_num = length(low_df)
high_num

df_sur <- df
df_sur <- subset(df_sur, futime != "NA")
df_sur = df_sur[df_sur$Row.names %in% c(high_df,low_df),]
dim(df_sur)
head(df_sur)
df_sur = cbind(df_sur,0)
title = "Overall Survival"
table(df_sur$fustat)
df_sur$futime = as.numeric(df_sur$futime) %/% 30 + 1
head(df_sur)
df_sur<-df_sur[!df_sur$futime>120,]    ## all
df_sur[df_sur$Row.names %in% high_df, 6] = "high"
df_sur[df_sur$Row.names %in% low_df,6] = "low"
df_sur[df_sur$futime==20,]
quantile(df[,gene])

df_sur$fustat = as.numeric(as.vector(df_sur$fustat))  # time
df_sur$futime = as.numeric(as.vector(df_sur$futime))  # OS
colnames(df_sur)[6] = "CLASS"
df_sur$CLASS = factor(df_sur$CLASS,levels = c("low","high"))
library(survival)
library("survminer")
mod = Surv(df_sur$futime, df_sur$fustat)
mfit = survfit(mod~df_sur$CLASS)
sur = survdiff(mod~df_sur$CLASS)
p.val <- 1 - pchisq(sur$chisq, length(sur$n) - 1)
p.val = signif(p.val,2)
p.val
head(df_sur)
dim(clinic)
quantile(mfit, probs=c(0.25, 0.5, 0.75), conf.int=FALSE)   #25%, 50% and 75% survival rate
# pdf("Figure_2E.pdf", width = 5, height = 5)
pdf("Figure_2F.pdf", width = 5, height = 5)
par(mar = c(4,4,2,1))
color <- c("#2b83ba","#d7191c")
plot(mfit,col = color, lwd = c(2,1.3,1.3,2,1.3,1.3),mark.time=T,conf.int = T,lty = c(1,3,3,1,3,3))
results_coxph = summary(coxph(mod~df_sur$CLASS))$coefficients
results_coxph_class = sub(pattern="df_sur\\$CLASS", replacement="", rownames(results_coxph))
results_coxph_hr = signif(results_coxph[1,"exp(coef)"],2)
results_coxph_hr_p = signif(results_coxph[1,"Pr(>|z|)"],2)
title(main = title,cex.main = 1.1,font.main = 1,line = 0.4)
xlab = "Months"
title(xlab = xlab,ylab="Percent survival",cex.lab=1.1,line = 2.5)
legend("topright",c(paste("Low",gene,"TPM",sep=" "),
                    paste("High",gene,"TPM",sep=" ")),
       col = color,lty = 1 ,lwd = 2,bty = "n",xjust = 1,
       y.intersp = 0.8,x.intersp = 0.5)
text(x = max(df_sur$futime) * 1.04,y = 0.73,
     labels = paste("Logrank p=",p.val,
                    "\n HR(",results_coxph_class,")=",results_coxph_hr,
                    "\n p(HR)=",results_coxph_hr_p,
                    "\n n(high)=",high_num,"\nn(low)=",low_num,
                    sep=""),pos = 2)
garbage <- dev.off()
#######


head(df)
# df<-df[!df$Row.names == "TCGA-06-0125",]
# df$fustat = ifelse(df$fustat == "Dead", "1", "0")
library(survivalROC)
rocCol=c("#d7191c","#2b83ba", "#abdda4")
aucText=c()
rt <- df
rt$futime <- rt$futime/365
# rt <- df_sur
head(rt)
dim(rt)
pdf(file="Figure_2H.pdf",width=5,height=5)
# pdf(file="Supplementary_Figure_3A.pdf",width=5,height=5)
# pdf(file="Supplementary_Figure_3B.pdf",width=5,height=5)
par(mar = c(4,4,2,1))
marker = rt[,3]
table(rt$futime > 2)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,3], predict.time =5, method="KM")
# roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,3], predict.time =3, method="KM")   ## GBM
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 2, cex.main=1.1, cex.lab=1.1, cex.axis=1.1, font=1.1)
aucText=c(aucText,paste0("five year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
# aucText=c(aucText,paste0("three year"," (AUC=",sprintf("%.3f",roc$AUC),")"))   ## GBM
abline(0,1)

roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,3], predict.time =3, method="KM")
# roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,3], predict.time =2, method="KM")   ## GBM
# aucText=c(aucText,paste0("two year"," (AUC=",sprintf("%.3f",roc$AUC),")"))    ## GBM
aucText=c(aucText,paste0("three year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)

roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,3], predict.time =1, method="KM")
aucText=c(aucText,paste0("one year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)

legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
dev.off()


