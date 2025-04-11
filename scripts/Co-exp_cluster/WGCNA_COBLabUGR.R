library(WGCNA)
enableWGCNAThreads(nThreads = 16)

# read counts data
PL_gene_counts <- read.csv("Outputs/PL_gene_counts.csv")

# Remove the first column containing the names of the probes (no for the numerical calculation))
rownames(PL_gene_counts)<- PL_gene_counts$geneID
PL_gene_counts$geneID <- NULL
datExpr0 = as.data.frame(t(PL_gene_counts));
names(datExpr0) = rownames(PL_gene_counts);
rownames(datExpr0) = names(PL_gene_counts);
rownames (datExpr0) = c("1", "2","3","4","5","6","7","8","9","10","11")

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

#  SAMPLE CLUSTERING ----
sampleTree = hclust(dist(datExpr0), method = "average");

plot(sampleTree, main = "Sample clustering to detect outliers in Placenta Tissue", sub="", xlab="",cex.lab = 1.5,cex.axis = 1.5, cex.main = 1.5, ) # sub & xlab to remove x-axis label

datExpr = as.data.frame(datExpr0)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

allTraits = PL_meta[, c(2, 25:26,41:45,47:49) ];

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
par(mfrow = c(1,2));
cex1 = 0.9;
# Adjusting plot axis
axis <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
i = 1
yaxis <- c(axis,1)
xaxis <- c(sft$fitIndices[,1],0)
# Scale-free topology fit index as a function of the soft-thresholding power
plot(xaxis, yaxis,
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# This line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 16;
adjacency = adjacency(datExpr, power = softPower);

#                     CALCULATE TOM-dissimilarity
#
#===============================================================================
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# PLOT TOM-dissimilarity

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

#                                Code chunk 6
#                   CLUSTER/MODULE COLORS AND GENE DENDROGRAM
# We like large modules, so we set the minimum module size relatively high
minModuleSize = 60;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

#                         PLOT DENDROGRAM WITH COLOURS
#                   REPRESENTING CO-EXPRESSION CLUSTERS/MODULES
#
#===============================================================================
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors in Pl Tissue")

#             MERGING MODULES WITH VERY SIMILAR PROFILES
#
#===============================================================================
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
plot(METree, main = "Clustering of module eigengenes in PL Tissue",
     xlab = "", sub = "")

#                          CHOOSING FUSION THRESHOLD
#
#===============================================================================
# Merging very similar co-expression clusters/modules
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules/clusters:
mergedMEs = merge$newMEs;

#                               COMPARATIVE OF
#                    DYNAMIC TREE CUT VS MERGED MODULES
#
#===============================================================================
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    cex.rowText = 1.5,
                    cex.colorLabels = 1, cex.dendroLabels = 0.9,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "FINE_02-dataInput.RData")

#                   QUANTIFY CLUSTER-TRAIT ASSOCIATIONS
#
#===============================================================================
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, allTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#    PLOT THE DATA OBTAINED FROM THE ASSOCIATIONS BETWEEN TRAITS AND CLUSTERS
#
#===============================================================================
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(7, 14, 3, 4));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(allTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               xLabelsAngle = 50,
               zlim = c(-1,1),
               cex.lab.x = 0.7,
               font.lab.x = 0.7,
               #              NO SIRVEN, COMPROBAR POR QUÉ:
               #               signifSymbols = c("***", "**", "*", ""),
               #               signifCutpoints = c(0, 0.001, 0.01, 0.05, 1),
               main = paste("Module-trait relationships in PL Tissue"))

#      RELATIONSHIP OF GENES WITH SIGNIFICANT CHARACTERISTICS AND CLUSTERS:
#              GENETIC SIGNIFICANCE (GS) AND MODULE MEMBERSHIP (MM)
#
#===============================================================================
# Define variable zscore containing the Z-score column of datTrait
allTraits$Zscore <- scale(allTraits$newborn_weight_kg)
zscore = as.data.frame(allTraits$Zscore);
names(zscore) = "zscore"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, zscore, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS_", names(zscore), sep="");
names(GSPvalue) = paste("p_GS_", names(zscore), sep="");

#                    IDENTIFYING GENES WITH HIGH GS AND MM
#
#===============================================================================
module = "purple"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for BW adjust by age and weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#         SUMMARY OUTPUT TO ANALYSE ONE OF THE CLUSTERS
#
#===============================================================================
# We have found clusters with high association with our trait of interest, and identified their core players through the Module Membership measure. We now merge this statistical information with the genetic annotation and write a file that summarizes the most important genetic annotation and results. It can be inspected in standard spreadsheet software such as Excel or Open Office Calc. Our expression data are only annotated by probe ID names


# To get the total number of genes present in our total data frame
names(datExpr)
names(datExpr)[moduleColors=="purple"]

# It informs us of the genes present in the "brown" cluter and which are present in the merged cluster.
# It will return probe IDs belonging to the brown cluster. To facilitate interpretation of the results, we used a probe annotation file provided by the manufacturer of the expression arrays to connect the probe IDs to universally recognized gene names and Entrez codes.

#          CREATE THE VARIABLE WITH THE NECESSARY INFORMATION THE GENES
#
#===============================================================================
# We now create a data frame that contains the following information for all probes:
# PROBE ID
# GENE SYMBOL
# LOCUS LINK ID (Entrez code),
# CLUSTER NAME
# BMI Z-SCORE AND GENE RELATIONSHIP
# MODULE MEMBERSHIP
# GENETIC SIGNIFICANCE

# Cluester will be ordered by their importance to the BMI Z-score, with the most significant on the left.

genes = names(datExpr)
# Create the starting data frame
geneInfo0 = data.frame(Code = genes,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order clusters by their significance for weight
modOrder = order(-abs(cor(MEs, zscore, use = "p")));
# Add cluster membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM_", modNames[modOrder[mod]], sep=""),
                       paste("p_MM_", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS_zscore));
geneInfo = geneInfo0[geneOrder, ]
geneInfo <- geneInfo0

rownames(geneInfo) <- gsub("\\/.*","",rownames(geneInfo))

write.csv(geneInfo, file = "geneInfo_PL.csv")

#  MODULE GENE ID  (GENES WHOSE MM > 0.8 WILL BE CONSIDERED HUBS)
#===============================================================================
table(gene0[,2] == geneInfo[,2])
COLOR <- names(table(geneInfo$moduleColor))
nombres <- colnames(geneInfo)
nombres <- gsub("MM_","",nombres)
nombres <- gsub("p_","",nombres)

for (i in 1:length(COLOR)) {
  color <- COLOR[i]
  tosave <- geneInfo[which(geneInfo$moduleColor == color),c(1:3,which(color == nombres))]
  a <- paste("MM_",color)
  a <- gsub(" ","",a)
  MM <- which(colnames(tosave)==a)
  tosave <- tosave[order(-tosave[,MM]),]
  colorname <- tosave[1,3]
  
  # SAVE DATA
  write.csv2(tosave,paste("EDA_PL/",colorname,".csv",sep=""))
}
