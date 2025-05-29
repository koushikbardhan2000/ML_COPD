# Load required library
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# Transpose expression data: samples as rows, genes as columns
datExpr <- t(batch_corrected)  # expr_matrix: genes × samples

# STEP 1: Check for good samples and genes
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# STEP 2: Sample clustering to check for outliers
sampleTree <- hclust(dist(datExpr), method = "average")
png("output/sample_clustering.png",  width = 6000, height = 2500, res = 300)
plot(sampleTree, main = "Sample Clustering to Detect Outliers", sub = "", xlab = "", cex.lab = 1, cex = 0.5)
dev.off()


# OPTIONAL: remove outliers if needed
cutHeight <- 100
clust <- cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 10)
keepSamples <- (clust == 1)
datExpr_clean <- datExpr[keepSamples, ]
pheno_hnVScopd_clean <- pheno_hnVScopd[keepSamples, ]

# STEP 3: Prepare phenotype data
# Ensure rownames match datExpr rows
traitData <- pheno_hnVScopd_clean
rownames(traitData) <- traitData$GSM_IDs
traitData <- traitData[rownames(datExpr), , drop = FALSE]
datTraits <- as.data.frame(traitData$phenotype)
colnames(datTraits) <- "Condition"

# STEP 4: Run blockwiseModules
net <- blockwiseModules(datExpr,
                        power = 8,
                        TOMType = "unsigned",
                        minModuleSize = 30,
                        reassignThreshold = 0,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "TOM",
                        verbose = 3)

# STEP 5: Convert numeric labels to colors
moduleColors <- labels2colors(net$colors)

# Plot gene dendrogram and module colors
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module Colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# STEP 6: Relate modules to traits (phenotype)
MEs <- net$MEs
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

# STEP 7: Plot heatmap of module-trait relationships
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = "Module-Trait Relationships")
