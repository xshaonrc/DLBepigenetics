library(data.table)
library(WGCNA)
library(flashClust)
setwd("./DLBproject/WGCNA")

commandIn <- commandArgs(F)
print(commandIn)

options(expressions = 500000)
inFilename <- "betas_1.csv" # DNA methylation beta value matrix

traits <- "Sample_Group,Gender,Age"

power.threshold <- 5 # predefined threshold 

cutline <- 400000

sd.top.quantile.cutoff <-
   0.95
# over 80% quantile; i.e. top20% variable CpGs
frac.subsample <- 1

minFraction.threshold <- 1 / 2

types <- "SampleBased"

topThres <- (1 - sd.top.quantile.cutoff) * 100

print(paste("power.threshold:", power.threshold))
print(paste("cutline:", cutline))
print(paste("sd.top.quantile.cutoff:", sd.top.quantile.cutoff))
print(paste("frac.subsample:", frac.subsample))
print(paste("minFraction.threshold:", minFraction.threshold))

methData <- fread(inFilename, data.table = F)
methbetaVals <- as.matrix(methData[,-c(1:5)])
rownames(methbetaVals) <- methData$ID
meta_filename <- "./data/SampleSheet.csv"

methData.matrix.traits <- fread(meta_filename, data.table = F)
rownames(methData.matrix.traits) <-
   methData.matrix.traits [, 1] ## the first column has to be the sample name;

traits.array <-
   unlist(strsplit(traits, split = ",")) ## have to be "," splited.
methData.matrix.selectedtraits <-
   methData.matrix.traits[, traits.array]

Groups <- rep(0,nrow(methData.matrix.selectedtraits))
Groups [methData.matrix.selectedtraits$Sample_Group == "DLB"] <- 1
methData.matrix.selectedtraits$Sample_Group <- Groups

###
rownames(methData.matrix.selectedtraits) <-
   methData.matrix.traits [, 1]

keepSamples <-
   colnames(methbetaVals) %in% rownames(methData.matrix.selectedtraits)
print(paste("Keeping samples with metadata:", sum(keepSamples)))
methData.matrix <- methbetaVals

# Filter based on minimum variance cutoff
if (sd.top.quantile.cutoff > 0) {
   # Filter for high variance sites
   methData.matrix.sd <-
      apply(methData.matrix[, keepSamples], 1, sd, na.rm = T)
   index <-
      methData.matrix.sd >= quantile(methData.matrix.sd,
                                     sd.top.quantile.cutoff,
                                     na.rm = T)
   methDat <- as.data.frame(t(methData.matrix[index, keepSamples]))
   colnames(methDat) <-
      rownames(methData.matrix[index, keepSamples])
   rownames(methDat) <-
      colnames(methData.matrix[index, keepSamples])
   
} else {
   print ("No sd quantile cutoff")
   methDat <-
      as.data.frame(t(methData.matrix[, keepSamples]))
   ## transpose: row is the sample ids, and the column is the CpG ids,
   colnames(methDat) <-  rownames(methData.matrix)
   rownames(methDat) <- colnames(methData.matrix[, keepSamples])
}

# Filter for minimum expression in dataset
# at least 50% of samples have some expression
keepCpGs <-
   (apply(methDat, 2, function (x)
      sum(!is.na(x)) / length(x)) == 1) 
print(paste(
   "Keeping CpGs with value in > 90 % of samples : ",
   sum(keepCpGs),
   "out of",
   length(keepCpGs)
))
methDat <- methDat[, keepCpGs]
# Subsample
# Matrix is transposed so samples are rows and CpG sites are cols
if (frac.subsample < 1.0) {
   print(paste("Subsample Fraction:", frac.subsample))
   selected <- sample(ncol(methDat), frac.subsample * ncol(methDat))
   write.table(selected, "wgcna.selectedcolumns.txt")
   methDat <- methDat[, selected]
} else {
   print("No subsampling - using all data")
}
widthscale <- max(nrow(methDat) / 100, 1.0)

gsg <-
   goodSamplesGenes(methDat, minFraction = minFraction.threshold, verbose = 3)
## more than half elements are not empty.
gsg$allOK

if (!gsg$allOK) {
   # Optionally, print the gene and sample names that were removed:
   if (sum(!gsg$goodGenes) > 0)
      printFlush(paste("Removing CpGs:",
                       paste(names(methDat)[!gsg$goodGenes], collapse = ", ")))
   
   if (sum(!gsg$goodSamples) > 0)
      printFlush(paste("Removing samples:",
                       paste(rownames(methDat)[!gsg$goodSamples], collapse = ", ")))
   
   # Remove the offending genes and samples from the data:
   methDat <- methDat[gsg$goodSamples, gsg$goodGenes]
}

sampleTree <- flashClust(dist(methDat), method = "average")


pdf(
   file = paste0(types, ".top", topThres, ".wgcna.Clustering.pdf"),
   width = 12 * widthscale,
   height = 9
)

par(cex = 0.6)

par(mar = c(0, 4, 2, 0))
plot(
   sampleTree,
   main = "Sample clustering to detect outliers",
   sub = "",
   xlab = "",
   cex.lab = 1.5,
   cex.axis = 1.5,
   cex.main = 2
)

# Cut the tree using cutline to eliminate outlier samples
# Plot a line to show the cut
abline(h = cutline, col = "red")

dev.off()


# Determine cluster under the line
clust <- cutreeStatic(sampleTree, cutHeight = cutline, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples <- (clust == 1)

print(paste("Samples Kept after Cutline: ", sum(keepSamples)))

datmeth <- data.matrix(methDat[keepSamples, ])
storage.mode(datmeth) <- "double"
nGenes <- ncol(datmeth)
nSamples <- nrow(datmeth)

datTraits <-
   data.matrix(methData.matrix.selectedtraits[match(rownames(datmeth), rownames(methData.matrix.selectedtraits)),]) 
   ## cannot handle string variable, currently; as string variable cannot be used for downstream correlation analysis.
storage.mode(datTraits) <- "double"
#traitColors=numbers2colors(as.numeric(datTraits[,1]), signed = FALSE)
traitColors <- numbers2colors((datTraits), signed = FALSE)

# Re-cluster samples
sampleTree2 <- flashClust(dist(datmeth), method = "average")

# Plot the sample dendrogram and the colors underneath.
pdf(
   file = paste0(types, ".top", topThres, ".wgcna.Clustering.trait.pdf"),
   width = 18 * widthscale,
   height = 24
)

par(cex = 0.6)

par(mar = c(0, 4, 2, 0))
plotDendroAndColors(sampleTree2,
                    traitColors,
                    groupLabels = colnames(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

### Automatic One-Step Network Construction and Module detection
# http://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-auto.pdf

enableWGCNAThreads()rownames(methData.matrix) <- CpGids;
colnames(methData.matrix) <-as.matrix(meth.SampleIds) ;

# Choose a set of soft-thresholding powers
powers <-
   c(seq(from = 1, to = 3, by = 0.5), c(4:10), seq(from = 12, to = 40, by =
                                                      2))
# Call the network topology analysis function
sft <-
   pickSoftThreshold(datmeth, powerVector = powers, verbose = 5)

power.threshold <- sft$powerEstimate

# Plot the results:
#sizeGrWindow(9, 5)
pdf(
   file = paste0(types, ".top", topThres, ".wgcna.PowerThreshold.pdf"),
   width = 12,
   height = 9
)

par(mfrow = c(1, 2))

cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(
   sft$fitIndices[, 1],
   -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
   xlab = "Soft Threshold (power)",
   ylab = "Scale Free Topology Model Fit,signed R^2",
   type = "n",
   main = paste("Scale independence")
)

text(
   sft$fitIndices[, 1],
   -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
   labels = powers,
   cex = cex1,
   col = "red"
)

# this line corresponds to using an R^2 cut-off of h
abline(h = 0.9, col = "black")
abline(v = power.threshold, col = "red")

# Mean connectivity as a function of the soft-thresholding power
plot(
   sft$fitIndices[, 1],
   sft$fitIndices[, 5],
   xlab = "Soft Threshold (power)",
   ylab = "Mean Connectivity",
   type = "n",
   main = paste("Mean connectivity")
)
text(
   sft$fitIndices[, 1],
   sft$fitIndices[, 5],
   labels = powers,
   cex = cex1,
   col = "red"
)
# this line corresponds to using an R^2 cut-off of h
abline(v = power.threshold, col = "red")

dev.off()


mergecutHeight.threshold <- 0.4
net <- blockwiseModules(
   datmeth,
   power = power.threshold,
   TOMType = "unsigned",
   minModuleSize = 30,
   reassignThreshold = 0,
   mergeCutHeight = mergecutHeight.threshold,
   numericLabels = TRUE,
   pamRespectsDendro = FALSE,
   saveTOMs = TRUE,
   saveTOMFileBase = paste0(types, ".TOM"),
   maxBlockSize = 40000,
   verbose = 3
)

# open a graphics window
#sizeGrWindow(12, 9)
pdf(
   file = paste0(
      types,
      ".top",
      topThres,
      ".power",
      power.threshold,
      ".mergeCutHeight",
      mergecutHeight.threshold,
      ".wgcna.ModuleDendro.pdf"
   ),
   width = 12 * widthscale,
   height = 9
)

# Convert labels to colors for plotting
mergedColors <- labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
   net$dendrograms[[1]],
   mergedColors[net$blockGenes[[1]]],
   "Module colors",
   dendroLabels = FALSE,
   hang = 0.03,
   addGuide = TRUE,
   guideHang = 0.05
)
dev.off()

## Save the modules and eigengene information
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs

geneTree <- net$dendrograms[[1]]

save(
   datmeth,
   datTraits,
   MEs,
   moduleLabels,
   moduleColors,
   geneTree,
   file = paste0(types, ".top5.wgcna.RData")
)

#### Relating modules to external information and identifying important genes
# http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-03-relateModsToExt.pdf

nSamples <- nrow(datmeth)

# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(datmeth, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

moduleTraitCor.table <-
   signif(cbind(moduleTraitCor, moduleTraitPvalue), 2)
colnames(moduleTraitCor.table) <-
   c(paste0(colnames(moduleTraitCor), ".cor"), paste0(colnames(moduleTraitPvalue), ".p"))
write.table(
   cbind(module = rownames(moduleTraitCor.table),
         moduleTraitCor.table),
   paste0(
      types,
      ".top",
      topThres,
      ".power",
      power.threshold,
      ".mergeCutHeight",
      mergecutHeight.threshold,
      ".wgcna.ModuleTraitCor.txt"
   ),
   quote = F,
   sep = "\t",
   row.names = F
)

pdf(
   file = paste0(
      types,
      ".top",
      topThres,
      ".power",
      power.threshold,
      ".mergeCutHeight",
      mergecutHeight.threshold,
      ".wgcna.ModuleTraitCor.pdf"
   ),
   width = 12,
   height = 9
)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2),
                   "\n(",
                   signif(moduleTraitPvalue, 1),
                   ")",
                   sep = "")

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

# Display the correlation values within a heatmap plot
labeledHeatmap(
   Matrix = moduleTraitCor,
   xLabels = colnames(datTraits),
   yLabels = names(MEs),
   ySymbols = names(MEs),
   colorLabels = FALSE,
   #colors = greenWhiteRed(50),
   colors = blueWhiteRed(50),
   textMatrix = textMatrix,
   setStdMargins = FALSE,
   cex.text = 0.5,
   zlim = c(-1, 1),
   main = paste("Module-trait relationships")
)

dev.off()


# Recalculate module eigengenes
MEs <- moduleEigengenes(datmeth, moduleColors)$eigengenes
MET <- orderMEs(cbind(MEs, datTraits))

# Plot the relationships among the eigengene modules and the trait
pdf(
   file = paste0(
      types,
      ".top",
      topThres,
      ".power",
      power.threshold,
      ".mergeCutHeight",
      mergecutHeight.threshold,
      ".wgcna.module.network.pdf"
   ),
   width = 12,
   height = 9
)

par(cex = 0.9)
plotEigengeneNetworks(
   MET,
   "",
   marDendro = c(0, 4, 1, 2),
   marHeatmap = c(3, 4, 1, 2),
   cex.lab = 0.8,
   xLabelsAngle = 90,
   plotAdjacency = F,
   printAdjacency = F
)
dev.off()
