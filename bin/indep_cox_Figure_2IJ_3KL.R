library(survival)

## survival_tcga_Figure_2EFH.R
head(df_sur)
df_sur <- df_sur[,c(4, 5, 3)]
head(df_cli)
df_cli <- df_cli[, c(4, 5, 6, 7, 8, 9, 10)]
colnames(df_cli)

sameSample=intersect(row.names(df_cli),row.names(df_sur))
df_sur=df_sur[sameSample,]
df_cli=df_cli[sameSample,]

rt=cbind(df_sur, df_cli)
head(rt)
# rt <- na.omit(rt)
rt$fustat <- as.numeric(rt$fustat)
rt$futime <- as.numeric(rt$futime)
rt <- subset(rt, age_at_initial_pathologic_diagnosis != "NA")
rt <- subset(rt, age_at_initial_pathologic_diagnosis != "")
# rt <- subset(rt, ldh1_mutation_found != "NA")
# rt <- subset(rt, ldh1_mutation_found != "")
rt <- subset(rt, radiation_therapy != "NA")
rt <- subset(rt, radiation_therapy != "")
rt <- subset(rt, neoplasm_histologic_grade != "")
rt$Age <- ifelse(rt$age_at_initial_pathologic_diagnosis <= 40,"<=40",">40")
rt$IDH <- rt$ldh1_mutation_found ## 逻辑
rt$'Radio' <- rt$radiation_therapy  ## 逻辑
rt$Gender <- rt$gender.demographic
head(rt)
colnames(rt)
table(rt$radiation_therapy)
table(rt$Radio)
table(rt$gender.demographic)
table(rt$neoplasm_histologic_grade)
rt$Grade <- rt$neoplasm_histologic_grade ## 逻辑
rt$Grade <- gsub("G[23]", "G2 G3", rt$neoplasm_histologic_grade)
table(rt$Grade)
colnames(rt)
rt_ind <- rt[, c(1, 2, 3, 11, 14, 15, 13)]   ## final tcga all
# rt_ind <- rt[, c(1, 2, 3, 11, 14, 13)]   ## final tcga sep
head(rt_ind)
table(rt$neoplasm_histologic_grade)
table(rt$disease_code)
table(rt_ind$Grade)
colnames(rt_ind)[7] <- "Treatment_Radio"    ## final tcga all
# colnames(rt_ind)[6] <- "Treatment_Radio"  ## final tcga sep
## data pre end

## cgga
# survival_cgga_Figure_3EFGHIJ.R
head(df)
colnames(df)
rt <- df
rt <- na.omit(rt)
rt$Grade <- gsub("^WHO III$", "Grade II III", rt$Grade)
rt$Grade <- gsub("^WHO II$", "Grade II III", rt$Grade)
rt$Grade <- gsub("^WHO IV$", "Grade IV", rt$Grade)
# rt <- rt[rt$Grade == "Grade IV", ]   ## sep GBM
rt <- rt[rt$Grade == "Grade II III", ]   ## sep LGG
rt$Age <- ifelse(rt$Age <= 40,"<=40",">40")
# rt$Age <- ifelse(rt$Age >= 41,">=41","<41")
rt$Gender  <- gsub("Female ", "Female", rt$Gender)
rt$Gender  <- gsub("Male ", "Male", rt$Gender)
rt$Gender  <- gsub("Male ", "Male", rt$Gender)
table(rt$`Radio_status (treated=1;un-treated=0)`)
rt$`Radio_status (treated=1;un-treated=0)`  <- ifelse(rt$`Radio_status (treated=1;un-treated=0)`, "Yes", "No")
table(rt$`Radio_status (treated=1;un-treated=0)`)
rt$`Chemo_status (TMZ treated=1;un-treated=0)` <- ifelse(rt$`Chemo_status (TMZ treated=1;un-treated=0)`, "Yes", "No")
table(rt$Grade)
table(rt$Age)
table(rt$Gender)
colnames(rt)
rt <- rt[, -1]
# rt_ind <- rt[, c(7, 8, 1, 6, 5, 4, 2, 11, 12, 13, 9, 10)]
rt_ind <- rt[, c(7, 8, 1, 6, 5, 2, 11, 12, 13, 9, 10)]  ## sep
colnames(rt_ind)
# colnames(rt_ind) <- c("futime", "fustat", "NID2", "Age", "Gender", "Grade", "Recurrent", "Mutaion_IDH", "Codeletion_1p19q", "Methylation_MGMTp", "Treatment_Radio", "Treatment_Chemo")
table(rt_ind$Treatment_Radio)
colnames(rt_ind) <- c("futime", "fustat", "NID2", "Age", "Gender", "Recurrent", "Mutaion_IDH", "Codeletion_1p19q", "Methylation_MGMTp", "Treatment_Radio", "Treatment_Chemo")   ## sep
## data pre end

## final
library("dplyr")
library(ggplot2)
library(survminer)
library(survival)
library(forestmodel)

panels <- list(
  list(width = 0.01),
  list(width = 0.1, display = ~variable, 
       fontface = "bold", heading = "Factor"),
  list(width = 0.1, display = ~level),
  list(width = 0.05, display = ~n, hjust = 1, heading = "Number"),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(width = 0.55, item = "forest", hjust = 0.5, 
       heading = "Hazard ratio", linetype = "dashed",line_x = 0),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(width = 0.12, display = ~ ifelse(reference, "Reference", sprintf("%0.2f (%0.2f, %0.2f)",
                                                                        trans(estimate), trans(conf.low), trans(conf.high))), display_na = NA),
  list(width = 0.05,
       display = ~ ifelse(reference, "", format.pval(p.value, digits = 1, eps = 0.001)),
       display_na = NA, hjust = 1, heading = "pvalue"),
  list(width = 0.03))

## multiForest
pdf("Figure_2J.pdf",width = 10, height = 5)
# pdf("Supplementary_Figure_3E.pdf",width = 10, height = 5)
# pdf("Supplementary_Figure_3F.pdf",width = 10, height = 5)
# pdf("Figure_3L.pdf",width = 10, height = 7)
# pdf("Supplementary_Figure_5E.pdf",width = 10, height = 7)
# pdf("Supplementary_Figure_5F.pdf",width = 10, height = 7)
# rt_ind <- rt_ind[, -7]
print(forest_model(coxph(Surv(futime, fustat) ~ ., rt_ind), factor_separate_line = T))
dev.off()

library(forestmodel)
colnames(rt_ind)[seq(3, ncol(rt_ind))]
covariates <- colnames(rt_ind)[seq(3, ncol(rt_ind))]
univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(futime, fustat)~', x)))
univ_formulas

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = rt_ind)})
pdf("Figure_2I.pdf",width = 10, height = 5)
# pdf("Supplementary_Figure_3C.pdf",width = 10, height = 5)
# pdf("Supplementary_Figure_3D.pdf",width = 10, height = 5)
# pdf("Figure_3K.pdf",width = 10, height = 7)
# pdf("Supplementary_Figure_5C.pdf",width = 10, height = 7)
# pdf("Supplementary_Figure_5D.pdf",width = 10, height = 7)
forest_model(model_list = univ_models, covariates = covariates, merge_models = T, panels)
## 6 x 8
dev.off()
df["CGGA_837", ]
