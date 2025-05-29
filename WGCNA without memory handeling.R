# Gene co-expression network analysis
# head(expr_hnVScopd,1)
expr_matrix <- expr_hnVScopd
pheno_data <- data.frame(
  SampleName = colnames(expr_hnVScopd),
  Condition = c(rep("HN", 161), rep("COPD", 92))  # Adjust based on your actual group sizes
)

# Load libraries
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()  # Enable multi-threading

# STEP 1: Transpose expression matrix (WGCNA needs samples as rows)
datExpr <- t(expr_matrix)

# Check for good samples and genes
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# STEP 2: Sample clustering to detect outliers
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

# STEP 3: Choose soft-thresholding power (β)
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot scale-free fit and mean connectivity
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n",
     main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=0.9, col="red")
abline(h=0.9, col="red")  # Usually pick power where R² > 0.9

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9, col="blue")

# STEP 4: Build network with selected power (e.g., 8 or as chosen)
softPower <- 8
adjacency <- adjacency(datExpr, power = softPower)

# STEP 5: Turn adjacency into topological overlap matrix (TOM)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# STEP 6: Cluster genes using TOM dissimilarity
geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)

# STEP 7: Identify modules with dynamic tree cut
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
table(dynamicMods)

# Convert numeric labels into colors
dynamicColors <- labels2colors(dynamicMods)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)

# STEP 8: Relate modules to phenotype
# Create design matrix from phenotype
traitData <- pheno_data
rownames(traitData) <- traitData$GSM  # Ensure rownames match colnames of expr_matrix
traitData <- traitData[colnames(expr_matrix), , drop = FALSE]  # Align samples
datTraits <- as.data.frame(traitData$Phenotype)
colnames(datTraits) <- "Condition"

# Module eigengenes
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

# Correlate with phenotype
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

# Heatmap of module-trait relationships
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
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
               main = "Module-phenotype relationships")
