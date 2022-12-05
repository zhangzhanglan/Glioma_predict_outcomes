library("limma")
library("WGCNA")

## 1F begin
#### clustering
sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = "Supplementary_Figure_1A.pdf", width = 14, height = 7)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, cex = 0.6)
####
abline(h = 133, col = "red")
dev.off()

#### cut above
clust = cutreeStatic(sampleTree, cutHeight = 133, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr0 = datExpr0[keepSamples, ]
dim(datExpr0)

#### clinical data
rownames(datExpr0)
temp_ggplot <- temp_ggplot[temp_ggplot$Row.names %in% rownames(datExpr0), ]
traitData <- data.frame(Control = ifelse(temp_ggplot$group == "control", 1, 0),
                        Glioma = ifelse(temp_ggplot$group == "glioma", 1, 0))
dim(data)
head(traitData)
datTraits <- traitData

#### clustering + clinical
sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(traitData, signed = FALSE)
pdf(file="Supplementary_Figure_1B.pdf",width=14,height=7)
plotDendroAndColors(sampleTree2, traitColors,
                    dendroLabels = FALSE, 
                    # groupLabels = names(datTraits),
                    groupLabels = names(traitData),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

#### power
allowWGCNAThreads() 
powers = c(1:20)
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
par(mfrow = c(1,2))
cex1 = 0.9
####
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red") 
####
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#### trans
sft    #### best power
softPower =sft$powerEstimate
adjacency = adjacency(datExpr0, power = softPower)
softPower

#### TOM
TOM = TOMsimilarity(adjacency)
dissTOM = 1- TOM

#### gene clustering
geneTree = hclust(as.dist(dissTOM), method = "average");
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)


####
minModuleSize = 25
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

####
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.6
abline(h=MEDissThres, col = "red")
####
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
pdf(file="Supplementary_Figure_1C.pdf", width = 9, height = 6)
plotDendroAndColors(geneTree, mergedColors,"Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
moduleColors = mergedColors
table(moduleColors)
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

####
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
pdf(file="Supplementary_Figure_1D.pdf",width=5.5,height=5.5)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

####
for (mod in 1:nrow(table(moduleColors))){ 
	modules = names(table(moduleColors))[mod]
	probes = colnames(datExpr0)
	inModule = (moduleColors == modules)
	modGenes = probes[inModule]
	write.table(modGenes, file =paste0("9_",modules,"_genes.txt"),sep="\t",row.names=F,col.names=F,quote=F)
}
table(moduleColors)
modules = names(table(moduleColors))[c(2, 5)]
modules
table(datTraits$Control)

modules = c("grey60", "plum2")
probes = colnames(datExpr0)
# inModule = (moduleColors == modules)
inModule = (moduleColors %in% modules)
modGenes = probes[inModule]
sigExp=t(datExpr0[,modGenes])
sigExpOut=rbind(ID=colnames(sigExp),sigExp)
sigExp[1:4, 1:4]
sigExp["22795_at", ]
write.table(sigExpOut, file="modGenes.exp.txt", sep="\t", quote=F, col.names=F)
