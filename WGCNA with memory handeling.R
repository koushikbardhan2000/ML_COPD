# Load required libraries
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# Transpose expression data: samples as rows, genes as columns
datExpr <- t(batch_corrected)

# STEP 1: Check for good samples and genes
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# STEP 2: Sample clustering to detect outliers
sampleTree <- hclust(dist(datExpr), method = "average")
png("output/sample_clustering.png", width = 6000, height = 2500, res = 300)
plot(sampleTree, main = "Sample Clustering to Detect Outliers", sub = "", xlab = "", cex.lab = 1, cex = 0.5)
dev.off()

# STEP 2.1: Remove known outliers
outliers_to_remove <- c("GSM349915", "GSM252799", "GSM924946", "GSM196992")
keepSamples <- !(rownames(datExpr) %in% outliers_to_remove)
datExpr_clean <- datExpr[keepSamples, ]
pheno_hnVScopd_clean <- pheno_hnVScopd[keepSamples, ]

# STEP 2.2: Re-check sample clustering after outlier removal
sampleTree2 <- hclust(dist(datExpr_clean), method = "average")
png("output/sample_clustering_clean.png", width = 6000, height = 2500, res = 300)
plot(sampleTree2, main = "Sample Clustering after cleaning Outliers", sub = "", xlab = "", cex.lab = 1, cex = 0.5)
dev.off()

# STEP 3: Prepare trait/phenotype data
traitData <- pheno_hnVScopd_clean
rownames(traitData) <- traitData$GSM_IDs
traitData <- traitData[rownames(datExpr_clean), , drop = FALSE]

# Encode phenotypes as numeric factor
datTraits <- data.frame(Condition = as.numeric(as.factor(traitData$phenotype)))
rownames(datTraits) <- rownames(datExpr_clean)  # Ensure alignment

# STEP 4: Construct network and identify modules
net <- blockwiseModules(datExpr_clean,
                        power = 8,
                        TOMType = "signed",
                        minModuleSize = 30,
                        reassignThreshold = 0,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "TOM",
                        verbose = 3)

# STEP 5: Convert numeric module labels to colors
moduleColors <- labels2colors(net$colors)

# Plot dendrogram with module colors
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module Colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# STEP 6: Relate modules to phenotype traits
MEs <- net$MEs
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr_clean))

# STEP 7: Plot module-trait heatmap
textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")")
dim(textMatrix) <- dim(moduleTraitCor)

png("output/tmp_sample.png", width = 6000, height = 2500, res = 300)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1, 1),
               main = "Module-Trait Relationships")
dev.off()
# Optional: Save gene-to-module assignments
geneModuleMembership <- data.frame(Gene = colnames(datExpr_clean), Module = moduleColors)




# STEP 1: Identify significantly associated modules (e.g., p < 0.05)
sigModules <- names(which(apply(moduleTraitPvalue, 1, function(p) any(p < 0.05))))
cat("Significant modules:", sigModules, "\n")

# STEP 2: Get gene-module assignments and expression matrix
geneModuleDF <- data.frame(Gene = colnames(datExpr_clean), Module = moduleColors)
selectedGenes <- list()

# STEP 3: For each significant module, find genes with high MM (Module Membership)
for (mod in sigModules) {
  modGenes <- geneModuleDF$Gene[geneModuleDF$Module == substring(mod, 3)]  # Remove "ME" prefix
  moduleColumn <- match(mod, colnames(MEs))  # Module eigengene column
  geneMMs <- abs(cor(datExpr_clean[, modGenes], MEs[, moduleColumn], use = "p"))  # Module Membership
  names(geneMMs) <- modGenes
  
  # Filter genes with high module membership (e.g., cor > 0.7)
  topGenes <- names(geneMMs)[geneMMs > 0.7]
  cat("Module", mod, ":", length(topGenes), "hub genes\n")
  selectedGenes[[mod]] <- topGenes
}

# STEP 4: Combine and export to file
importantGenes <- unique(unlist(selectedGenes))
write.csv(importantGenes, file = "output/potential_biomarkers.csv", row.names = FALSE)

cat("Saved", length(importantGenes), "potential biomarkers to 'output/potential_biomarkers.csv'\n")
