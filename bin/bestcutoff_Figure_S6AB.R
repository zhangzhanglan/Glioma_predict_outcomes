
## survival_tcga_Figure_2EFGH.R 
head(df)
library(survival)
library(survminer)
df$fustat <- as.numeric(df$fustat)
df$futime <- as.numeric(df$futime)
sur.cut <- surv_cutpoint(df,   time= 'futime',
                         event = 'fustat' ,
                         variables = 'NID2' )
sur.cut
sur.cat <- surv_categorize(sur.cut)
table(sur.cat$NID2)

sur.cut.plot <- surv_cutpoint(
  df,
  time = "futime",
  event = "fustat",
  variables = "NID2")
summary(sur.cut.plot)

pdf("Supplementary_Figure_6A.pdf", width = 5,height = 5, onefile=FALSE)
plot(sur.cut.plot, "NID2", palette = c("#d7191c","#2b83ba"))
dev.off()

fit <- survfit(Surv(futime, fustat) ~NID2, data =sur.cat)
ggpar( ggsurvplot(
  fit,   data = sur.cat, pval = TRUE,           
  conf.int = TRUE)
)
sur.cut$cutpoint$cutpoint
## [1] 1.163174
## [1] 2.513645  log2(tpm + 1)
## [1] 2.266865  log2(tpm)
## [1] 0.7566827  log10(tpm + 1)  #drug

gene <- "NID2"
head(rt)
head(df)
high_df = rownames(rt[rt[,gene] > sur.cut$cutpoint$cutpoint ,1,drop = F])
low_df = rownames(rt[rt[,gene] <= sur.cut$cutpoint$cutpoint,1,drop = F])
high_num = length(high_df)
low_num = length(low_df)
high_num

df_sur <- df
head(df_sur)
df_sur$Row.names <- rownames(df_sur)
df_sur <- subset(df_sur, futime != "NA")
df_sur = df_sur[df_sur$Row.names %in% c(high_df,low_df),]
head(df_sur)
df_sur = cbind(df_sur,0)
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

pdf("Supplementary_Figure_6B.pdf", width = 5,height = 5)
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
