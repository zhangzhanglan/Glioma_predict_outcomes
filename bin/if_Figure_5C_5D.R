rbg <- exp <- read.csv("../IF/21-20492-5-1.csv")   ## final GBM
# rbg <- exp <- read.csv("../IF/20-05150-3.csv")   ## final LGG
head(rbg)
plot(x = rbg$R, y = rbg$G)
library(ggplot2)
library(dplyr)
library(ggpubr)
cor(rbg$R,rbg$G,method="pearson")
# pdf(file="Figure_5D_up.pdf", width=3, height=3)
pdf(file="Figure_5D_do.pdf", width=3, height=3)
ggplot(rbg, aes(x=R, y=G))+ geom_point(size=0.1, shape=15)+geom_smooth(method=lm) + xlab("NID2 intensity (A.U.)") + ylab("CD31 intensity (A.U.)") + theme_bw()+ stat_cor(method = 'spearman', aes(x =R, y =G))
dev.off()

library(reshape)
head(rbg)
df <- data.frame(x = seq_along(rbg[, 1]), rbg)
head(df)
df <- df[, -4]
df <- melt(df, id.vars = "x")
cols <- c("#D43F3A", "#5CB85C")
# pdf(file="Figure_5C_up.pdf", width=4, height=3)
pdf(file="Figure_5C_do.pdf", width=4, height=3)
ggplot(df, aes(x = x, y = value, color = variable)) +
  geom_line() +
  scale_color_manual(values = cols, labels =c("NID2", "CD31")) + 
  xlab("Profile (nm)") + 
  ylab("Fluorescence intensity (A.U.)") +
  theme_bw() +
  guides(color = guide_legend(title = "Channel")) 
dev.off()
