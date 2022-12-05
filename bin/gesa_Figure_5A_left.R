library("RColorBrewer")
# install.packages("GGally")
# install.packages("scales")
# devtools::install_github("briatte/ggnet")
library(GGally)
library(ggnet)
library(network)
library(sna)
library(ggplot2)
library(igraph)
library(knitr)
library(kableExtra)
bip = network(h.dat,
              matrix.type = "bipartite",
              ignore.eval = FALSE,
              vertex.attrnames = c("Hallmarks", "Genes"),
              names.eval = "weights")
x = network.vertex.names(bip)
x
type <- c(rep("Genes", times=116), rep("Hallmarks", times = 20))

## start plot
library(dplyr)
library(tidyr)
library(ggplot2) # needs to be version ≥ 2.1.0
library(scales)
library(geomnet)
# install.packages('geomnet')
# devtools::install_github("sctyner/geomnet")
x = network.vertex.names(bip)
# x[117:136] <- seq(1, 20, by = 1)  ## hall is no.
hall <- x[117:136]  ## hall is pathway
hall
head(hall)
length(x)
rep("white",116)
bip %v% "Module" = c(rep("Genes", 116), "M5:Transmenbrane transporter", "M3:Vessel morphogensis", "M1:Cell cycle regulation", "M2:Immune response", "M2:Immune response", "M2:Immune response", "M2:Immune response", "M1:Cell cycle regulation", "M1:Cell cycle regulation", "M2:Immune response", "M2:Immune response", "M2:Immune response", "M2:Immune response", "M5:Transmenbrane transporter", "M5:Transmenbrane transporter", "M4:External encapsulating", "M5:Transmenbrane transporter", "M1:Cell cycle regulation", "M5:Transmenbrane transporter", "M5:Transmenbrane transporter")
bip %v% "Hallmarks" = type
type
set.seed(576)
RColorBrewer::brewer.pal(8, "Paired")[2]
colors_6<-c("#3cb346","#00abf0","#d75427",
            "#2e409a","#942d8d","#eeb401")
mm.col <- c("Genes" = "#A6CEE3", "M1:Cell cycle regulation" = "#3cb346", "M2:Immune response" = "#d75427",
            "M3:Vessel morphogensis" = "#1F78B4", "M4:External encapsulating" = "#942d8d", 
            "M5:Transmenbrane transporter" = "#eeb401")

ggnet.cichlids <- 
  ggnet2(bip, 
         node.size = 12, 
         node.color = "Module", 
         color.palette = mm.col,
         color.legend = "Module", 
         edge.size = 1, 
         edge.color = "grey50",
         arrow.size = 9,  arrow.gap = 0.027, 
         legend.size=20, label = x, label.size = 2)  +
  guides(color=guide_legend(keyheight=0.5, default.unit="inch",
                            override.aes = list(size=6)))
pdf("Figure_5A_l.pdf", width=14, height=7)
ggnet.cichlids
dev.off()

