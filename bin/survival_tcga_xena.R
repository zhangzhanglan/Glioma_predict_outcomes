library(data.table)
library(dplyr)
library(tidyr)
# phenotype <- fread('TCGA-LGG.GDC_phenotype.tsv.gz', header=T)
# phenotype <- fread('TCGA-GBM.GDC_phenotype.tsv.gz', header=T)
phenotype<-XenaGenerate(subset = XenaCohorts =="GDC TCGA Glioblastoma (LGG)")%>% 
  XenaFilter(filterDatasets    = "TCGA-LGG.GDC_phenotype.tsv") %>% 
  XenaQuery() %>%
  XenaDownload() %>% 
  XenaPrepare() 
phenotype<-XenaGenerate(subset = XenaCohorts =="GDC TCGA Glioblastoma (GBM)")%>% 
  XenaFilter(filterDatasets    = "TCGA-GBM.GDC_phenotype.tsv") %>% 
  XenaQuery() %>%
  XenaDownload() %>% 
  XenaPrepare() 

# 提取肿瘤的表型信息
tumor.sample <- sapply(as.character(phenotype$submitter_id.samples), function(x){
  number <- as.numeric(unlist(strsplit(unlist(strsplit(as.character(x), split="-"))[4], split = "[A-Z]"))[1])
  if(number<=9){
  # if(number<=11){
    id <- as.character(x)
    return(id)
  }
})
tumor.sample.id <- unlist(tumor.sample)
tumor.phenotype <- phenotype[match(tumor.sample.id, phenotype$submitter_id.samples), ]
tumor.phenotype <- as.data.frame(tumor.phenotype)
rownames(tumor.phenotype) <- tumor.phenotype$submitter_id.samples
# normal <- c("TCGA-06-0681-11A", "TCGA-06-AABW-11A", "TCGA-06-0678-11A", "TCGA-06-0675-11A", "TCGA-06-0680-11A") ## GBM
# table(rownames(tumor.phenotype) %in% normal)
table(tumor.phenotype$radiation_therapy)
table(tumor.phenotype$gender.demographic)
table(tumor.phenotype$age_at_initial_pathologic_diagnosis)
table(tumor.phenotype$inherited_genetic_syndrome_result)    ## 1p19q codeleted 
table(tumor.phenotype$ldh1_mutation_found)
table(tumor.phenotype$neoplasm_histologic_grade)
table(tumor.phenotype$primary_diagnosis.diagnoses)
table(tumor.phenotype$progression_or_recurrence.diagnoses)
table(tumor.phenotype$disease_code)
# query <- c("age_at_initial_pathologic_diagnosis", "gender.demographic", "radiation_therapy", "inherited_genetic_syndrome_result", "ldh1_mutation_found", "neoplasm_histologic_grade", "primary_diagnosis.diagnoses", "disease_code")     ## LGG
query <- c("age_at_initial_pathologic_diagnosis", "gender.demographic", "radiation_therapy", "ldh1_mutation_found", "neoplasm_histologic_grade", "primary_diagnosis.diagnoses", "disease_code")     ## LGG
query <- c("age_at_initial_pathologic_diagnosis", "gender.demographic", "radiation_therapy", "primary_diagnosis.diagnoses", "disease_code")

tumor.phenotype[1:2, 1:2]
uniq.sample.id <- unique(tumor.phenotype$submitter_id)
uniq.tumor.phenotype <- tumor.phenotype[match(uniq.sample.id, tumor.phenotype$submitter_id), ]
rownames(uniq.tumor.phenotype) <- uniq.tumor.phenotype$submitter_id
uniq.tumor.phenotype[1:2, 1:2]


## surivival
# uniq.tumor.phenotype$year_of_death.demographic
table(uniq.tumor.phenotype$vital_status.demographic)
table(uniq.tumor.phenotype$submitter_id.samples %in% normal)
uniq.tumor.phenotype <- subset(uniq.tumor.phenotype, vital_status.demographic != "Not Reported")
uniq.tumor.phenotype <- subset(uniq.tumor.phenotype, vital_status.demographic != "")
uniq.tumor.phenotype$OS<- ifelse((uniq.tumor.phenotype$vital_status.demographic == "Alive"), uniq.tumor.phenotype$days_to_last_follow_up.diagnoses,uniq.tumor.phenotype$days_to_death.demographic)
uniq.tumor.phenotype$OS_status<- as.factor(uniq.tumor.phenotype$vital_status.demographic)
uniq.tumor.phenotype <- select(uniq.tumor.phenotype, "OS_status", "OS")
head(uniq.tumor.phenotype)

## clinical phenotype
uniq.tumor.phenotype <- dplyr::select(uniq.tumor.phenotype, query)
# phe <- uniq.tumor.phenotype
lgg_phe <- uniq.tumor.phenotype
gbm_phe <- uniq.tumor.phenotype
colnames(gbm_phe)
colnames(lgg_phe)
table(lgg_phe$ldh1_mutation_found)
table(lgg_phe$neoplasm_histologic_grade)
dim(gbm_phe)
gbm_phe$ldh1_mutation_found <- rep("NA", 606)
gbm_phe$neoplasm_histologic_grade <- rep("G4", 606)

phe <- rbind(lgg_phe, gbm_phe)
table(phe$ldh1_mutation_found)
colnames(phe)
rownames(phe)
dim(phe)
save(phe, file = "../data/survival_tcga_xena.Rdata")
# save(phe, file = "../data/clinical_tcga_xena.Rdata")
# phe <- gbm_phe

table(lgg_phe$gender.demographic)
table(gbm_phe$gender.demographic)
lgg_phe$Age <- ifelse(lgg_phe$age_at_initial_pathologic_diagnosis <= 40,"<=40",">40")
table(lgg_phe$Age)
gbm_phe$Age <- ifelse(gbm_phe$age_at_initial_pathologic_diagnosis <= 40,"<=40",">40")
table(gbm_phe$Age)

table(lgg_phe$radiation_therapy)
table(gbm_phe$radiation_therapy)

table(lgg_phe$ldh1_mutation_found)
table(gbm_phe$ldh1_mutation_found)

table(lgg_phe$neoplasm_histologic_grade)
table(gbm_phe$neoplasm_histologic_grade)

table(lgg_phe$primary_diagnosis.diagnoses)
table(gbm_phe$primary_diagnosis.diagnoses)

table(phe$gender.demographic)
phe$Age <- ifelse(phe$age_at_initial_pathologic_diagnosis <= 40,"<=40",">40")
table(phe$Age)

table(phe$radiation_therapy)

table(phe$ldh1_mutation_found)

table(phe$neoplasm_histologic_grade)

table(phe$primary_diagnosis.diagnoses)
