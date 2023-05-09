## survival_tcga_Figure_2EFH.R
data[1:4, 1:4]
dim(data)
ins <- c("VEGFA", "FLT4", "KDR", "NRP1", "VWF", "PECAM1", 'ENG',"DLL4", "NID2")
exp <-data[row.names(data) %in% ins, ]
rownames(exp)


## angio
library(tidyr)
library(dplyr)
library(tibble)
rownames(exp)
exp[1:4, 1:4]
exp <- exp[sort(rownames(exp)), ]
rownames(exp)
rownames(exp) <- c("DLL4", "CD105 (ENG)", "VEGFR3 (FLT4)", "VEGFR2 (KDR)", "NID2", "NRP1", "CD31 (PECAM1)", "VEGFA", "von Willebrand factor (VWF)")
exp_temp <- as.data.frame(t(exp))
exp_longer<- exp_temp %>% 
  rownames_to_column(var = 'Sample') %>% 
  pivot_longer( cols =   c("DLL4", "CD105 (ENG)", "VEGFR3 (FLT4)", "VEGFR2 (KDR)", "NRP1", "CD31 (PECAM1)", "VEGFA", "von Willebrand factor (VWF)"),
                names_to = 'Gene',
                values_to = 'Expression')
head(exp_longer)

set.seed(123)

library(dplyr)
library(ggplot2)
library(ggpubr)
p <- ggplot(exp_longer, aes(x = NID2, y = Expression)) +
  geom_point() + 
  geom_smooth(method = "lm", col = "blue") +
  facet_wrap(~ Gene, nrow = 2) +
  stat_cor(method = "pearson", label.y = 9.5) +
  ylab("The expression of vascular-related marker genes") +
  xlab("NID2 expression")
pdf("Figure_5E.pdf",height=5, width=14)
print(p)
dev.off()
