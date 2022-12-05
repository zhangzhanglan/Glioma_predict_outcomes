library(tidyverse)

stad_fpkm <- read_tsv("TCGA-LGG.htseq_fpkm.tsv.gz")
# stad_fpkm <- read_tsv("TCGA-GBM.htseq_fpkm.tsv.gz")

gtf_v22 <- read_tsv("gencode.v22.annotation.gene.probeMap") %>% dplyr::select(1,2)
colnames(stad_fpkm)[1] <- "id"
stad_fpkm1 <- inner_join(gtf_v22, stad_fpkm, by = "id") %>% dplyr::select(-1)
stad_fpkm1[1:4, 1:4]
stad_fpkm1 <- aggregate(.~ gene, stad_fpkm1, mean)
stad_fpkm1 <- stad_fpkm1 %>% column_to_rownames("gene")
stad_fpkm1 <- 2^stad_fpkm1-1

fpkmToTpm <- function(fpkm){ exp(log(fpkm) - log(sum(fpkm)) + log(1e6))}
stad_tpm <- apply(stad_fpkm1, 2, fpkmToTpm)

save(stad_tpm, file = "TCGA-LGG.htseq_fpkmtotpm.Rdata")
# save(stad_tpm, file = "TCGA-GBM.htseq_fpkmtotpm.Rdata")
## end

colnames(stad_tpm)
# stad_tpm <- log2(stad_tpm + 0.001)
# stad_tpm <- log2(stad_tpm + 1)
# stad_tpm <- log10(stad_tpm + 1)