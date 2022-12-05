library(sva)
library(limma)
## http://www.cgga.org.cn/download.jsp
file1="CGGA.mRNAseq_693.RSEM-genes.20200506.txt"  
file2="CGGA.mRNAseq_325.RSEM-genes.20200506.txt"

rt1=read.table(file1,sep="\t",header=T,check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
data1=matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames1)
data1=log2(data1+1)

rt2=read.table(file2,sep="\t",header=T,check.names=F)
rt2=as.matrix(rt2)
rownames(rt2)=rt2[,1]
exp2=rt2[,2:ncol(rt2)]
dimnames2=list(rownames(exp2),colnames(exp2))
data2=matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)
data2=log2(data2+1)

dim(data1)
dim(data2)

sameGene=intersect(row.names(data1),row.names(data2))
data=cbind(data1[sameGene,],data2[sameGene,])
dim(data)

batchType=c(rep(1,ncol(data1)),rep(2,ncol(data2)))
data=ComBat(data, batchType, par.prior=TRUE)

ins <- c("NID1", "NID2")
exp <-data[row.names(data) %in% ins, ]
exp[1:2, 1:2]
dim(exp)
exp["NID2", "CGGA_837"]

clinic1<-read.table("CGGA.mRNAseq_325_clinical.20200506.txt",header = T,check.names = F,sep = "\t")
# clinic3<-read.table("CGGA.Methyl_array_159_clinical.20200506.txt",header = T,check.names = F,sep = "\t")
# clinic4<-read.table("CGGA.WEseq_286_clinical.20200506.txt",header = T,check.names = F,sep = "\t")
clinic2<-read.table("CGGA.mRNAseq_693_clinical.20200506.txt",header = T,check.names = F,sep = "\t")
clinic <- rbind(clinic1, clinic2)
table(clinic$Grade)
# name <- intersect(clinic$CGGA_ID, clinic3$CGGA_ID)
# length(name)
# name <- intersect(name, clinic4$CGGA_ID)
# name <- intersect(clinic3$CGGA_ID, clinic4$CGGA_ID)
colnames(clinic)
# age <- clinic$Age
# age <- na.omit(age)
# mean(age)
table(clinic$IDH_mutation_status)
# dim(clinic4)
table(clinic$Grade)
#### clinic all pre
# clinic_II <- clinic[clinic$Grade=="WHO II" & clinic$IDH_mutation_status == "Mutant", c(1,7,8)]    ## G
# clinic_II <- clinic[clinic$Grade=="WHO II" & clinic$IDH_mutation_status == "Wildtype", c(1,7,8)]    ## H
# clinic_III <- clinic[clinic$Grade=="WHO III" & clinic$IDH_mutation_status == "Mutant", c(1,7,8)]   ## G
# clinic_III <- clinic[clinic$Grade=="WHO III" & clinic$IDH_mutation_status == "Wildtype", c(1,7,8)]    ## H
# clinic_II<-clinic[clinic$Grade=="WHO II", c(1,7,8)]
# clinic_III<-clinic[clinic$Grade=="WHO III", c(1,7,8)]
clinic_II<-clinic[clinic$Grade=="WHO II", ]
clinic_III<-clinic[clinic$Grade=="WHO III", ]
clinic<-rbind(clinic_II, clinic_III)  ## LGG  FGH
# clinic <- clinic[clinic$Grade=="WHO IV", c(1,7,8)]    ## GBM
clinic <- clinic[clinic$Grade=="WHO IV", ]    ## GBM
clinic <- clinic[, c(1,7,8)]    ## ALL  survival
# head(clinic)
clinic<-clinic[!is.na(clinic$OS),]
clinic<-clinic[!is.na(clinic$`Censor (alive=0; dead=1)`),]
# clinic<-clinic[!clinic$OS<30,]
clinic<-clinic[!is.na(clinic$CGGA_ID),]

rownames(clinic) <- clinic$CGGA_ID
clinic <- clinic[, -1]
head(clinic)
## clinic prepare

colnames(clinic) <- c("futime", "fustat")
rt <- exp
rt <- t(rt)
sample <- intersect(rownames(rt), rownames(clinic))
head(sample)
rt = rt[sample,]
clinic =clinic[sample,]
df <- merge(rt, clinic, by= "row.names", all.x=TRUE)
rownames(df) <- df$Row.names
df <- df[, -1]
head(df)
dim(df)
## df prepare

##### 
head(df)
highcutoff <- 50
lowcutoff <- 50
gene <- "NID2"
high_df = rownames(rt[rt[,gene] > quantile(rt[,gene],as.numeric(highcutoff)/100),1,drop = F])
low_df = rownames(rt[rt[,gene] < quantile(rt[,gene],as.numeric(lowcutoff)/100),1,drop = F])
high_num = length(high_df)
low_num = length(low_df)
high_num

df_sur <- df
df_sur <- subset(df_sur, futime != "NA")
df_sur = df_sur[df_sur$Row.names %in% c(high_df,low_df),]
head(df_sur)
df_sur = cbind(df_sur, 0)
df_sur[, 4]
title = "Overall Survival"

df_sur$futime = as.numeric(df_sur$futime) %/% 30 + 1
df_sur[df_sur$Row.names %in% high_df, 6] = "high"
df_sur[df_sur$Row.names %in% low_df,6] = "low"
head(df_sur)

df_sur = as.data.frame(df_sur)
df_sur$fustat = as.numeric(as.vector(df_sur$fustat))  # time
df_sur$futime = as.numeric(as.vector(df_sur$futime))  # OS
colnames(df_sur)[6] = "CLASS"
df_sur$CLASS = factor(df_sur$CLASS,levels = c("low","high"))
head(df_sur)


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
pdf(file = "Figure_3E.pdf", width=5, height=5)
pdf(file = "Figure_3F.pdf", width=5, height=5)
pdf(file = "Figure_3G.pdf", width=5, height=5)
pdf(file = "Figure_3H.pdf", width=5, height=5)
pdf(file = "Figure_3I.pdf", width=5, height=5)
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
library(survivalROC)
rocCol=c("#d7191c","#2b83ba", "#abdda4")
aucText=c()
rt <- df
rt$futime <- rt$futime/365
head(rt)
dim(rt)
pdf(file = "Figure_3J.pdf", width=5, height=5)
pdf(file = "Supplementary_Figure_5A.pdf", width=5, height=5)
pdf(file = "Supplementary_Figure_5B.pdf", width=5, height=5)
par(mar = c(4,4,2,1))
marker = rt[,3]
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,3], predict.time =5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1],
     # main = "Time-dependent ROC curve", 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 2, cex.main=1.1, cex.lab=1.1, cex.axis=1.1, font=1.1)
aucText=c(aucText,paste0("five year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
abline(0,1)

roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,3], predict.time =3, method="KM")
aucText=c(aucText,paste0("three year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)

roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,3], predict.time =1, method="KM")
aucText=c(aucText,paste0("one year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)

legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
dev.off()

## table
colnames(clinic)
tbl <- prop.table(table(clinic$Grade))
cbind(tbl,prop.table(tbl))
table(tble$`Radio_status (treated=1;un-treated=0)`)

tblFun <- function(x){
  tbl <- table(x)
  res <- cbind(tbl,round(prop.table(tbl)*100,2))
  colnames(res) <- c('Count','Percentage')
  res
}
tble <- clinic
colnames(tble)
tble$Age <- ifelse(tble$Age <= 40,"<=40",">40")
tble$Gender  <- gsub("Female ", "Female", tble$Gender)
tble$Gender  <- gsub("Male ", "Male", tble$Gender)
tble$Gender  <- gsub("Male ", "Male", tble$Gender)
tble <- tble[, -7]
tble <- tble[, -7]
do.call(rbind,lapply(tble[2:11], tblFun))
