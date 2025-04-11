#               CREATION OF DATAFRAME WITH GENE EXPRESSION DATA
#                       AND ITS ADJUSTMENT FOR ANALYSIS
#
#===============================================================================
# We work with two sets:
nSets = 2;
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("PL tissue", "CB tissue")
shortLabels = c("PL", "CB")
# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)

multiExpr[[1]] = list(data = as.data.frame(t(PL_gene_counts[])));
names(multiExpr[[1]]$data)
rownames(multiExpr[[1]]$data)
multiExpr[[2]] = list(data = as.data.frame(t(CB_gene_counts[])));
names(multiExpr[[2]]$data)
rownames(multiExpr[[2]]$data)
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)
checkSets(multiExpr, checkStructure = T)

gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}
#                              SAMPLE CLUSTERING
#
#===============================================================================
sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);

# Check the size of the leftover data
exprSize = checkSets(multiExpr)
exprSize

# Form a multi-set structure that will hold the clinical traits
# Initialize the Traits list to hold the clinical data
Traits = vector(mode = "list", length = 2)  # Length should be 2, for PL and CB datasets

# Loop over each gene expression dataset (PL_gene_counts and CB_gene_counts)
for (set in 1:2) {
  # Get the sample names (rownames) from the current gene expression dataset
  setSamples = rownames(multiExpr[[set]]$data)
  
  # Remove the "PL" or "CB" prefix and convert to numeric IDs
  setSamplesNumeric = as.numeric(sub("^(PL|CB)", "", setSamples))
  
  # Match these numeric sample IDs with the study_ID in allTraits
  traitRows = match(setSamplesNumeric, allTraits$study_ID)
  
  # Check if any samples didn't match (e.g., NA values)
  if (any(is.na(traitRows))) {
    warning("Some samples in multiExpr[[", set, "]] don't have matching study_IDs in allTraits.")
    print(setSamples[is.na(traitRows)])  # Print unmatched sample names for debugging
  }
  
  # Assign the clinical data (excluding the study_ID column)
  Traits[[set]] = list(data = allTraits[traitRows, -1])  # Exclude the first column (study_ID)
  
  # Set the row names of the clinical data to be the corresponding study_ID values
  rownames(Traits[[set]]$data) = allTraits$study_ID[traitRows]
}

collectGarbage();
# Define data set dimensions
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;
nSets = checkSets(multiExpr)$nSets
save(multiExpr, Traits, nGenes, nSamples, setLabels, shortLabels, exprSize,
     file = "WGCNA_Consensus_Input.RData");

# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]]);
collectGarbage();
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
#                  SET ESTABLISH THE POWER OF THE SOFT THRESHOLD
#                 (based on the graph, where the red line cuts off)
softPower = c(12,16);
# Initialize an appropriate array to hold the adjacencies
adjacencies = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate adjacencies in each individual data set
for (set in 1:nSets)
  adjacencies[set, , ] = abs(cor(multiExpr[[set]]$data, use = "p"))^softPower[set];

#                     CALCULATE TOM-dissimilarity
#
#===============================================================================
# Initialize an appropriate array to hold the TOMs
TOM = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate TOMs in each individual data set
for (set in 1:nSets)
  TOM[set, , ] = TOMsimilarity(adjacencies[set, , ]);
# Define the reference percentile
scaleP = 0.95
# Set RNG seed for reproducibility of sampling | set.seed(12345)

# Sample sufficiently large number of TOM entries
nSamples = as.integer(1/(1-scaleP) * 1000);
# Choose the sampled TOM entries
scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
TOMScalingSamples = list();
# These are TOM values at reference percentile
scaleQuant = rep(1, nSets)
# Scaling powers to equalize reference TOM values
scalePowers = rep(1, nSets)
# Loop over sets
for (set in 1:nSets)
{
  # Select the sampled TOM entries
  TOMScalingSamples[[set]] = as.dist(TOM[set, , ])[scaleSample]
  # Calculate the 95th percentile
  scaleQuant[set] = quantile(TOMScalingSamples[[set]],
                             probs = scaleP, type = 8);
  # Scale the male TOM
  if (set>1)
  {
    scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
    TOM[set, ,] = TOM[set, ,]^scalePowers[set];
  }
}
#                           PLOT TOM-dissimilarity
#
#===============================================================================
# For plotting, also scale the sampled TOM entries
scaledTOMSamples = list();
for (set in 1:nSets)
  scaledTOMSamples[[set]] = TOMScalingSamples[[set]]^scalePowers[set]
# Open a suitably sized graphics window
sizeGrWindow(6,6)
#pdf(file = "Plots/TOMScaling-QQPlot.pdf", wi = 6, he = 6);
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[2]], plot.it = TRUE, cex = 0.6,
                    xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[2]),
                    main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[2]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20);
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
#dev.off();
consensusTOM = pmin(TOM[1, , ], TOM[2, , ]);

# Clustering
consTree = hclust(as.dist(1-consensusTOM), method = "average");
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 60;
# Module identification using dynamic tree cut:
unmergedLabels = cutreeDynamic(dendro = consTree, distM = 1-consensusTOM,
                               deepSplit = 2, cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamRespectsDendro = FALSE );
unmergedColors = labels2colors(unmergedLabels)

plotDendroAndColors(consTree, unmergedColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#             MERGING MODULES WITH VERY SIMILAR PROFILES
#
#===============================================================================
# Calculate module eigengenes
unmergedMEs = multiSetMEs(multiExpr, colors = NULL, universalColors = unmergedColors)
# Calculate consensus dissimilarity of consensus module eigengenes
consMEDiss = consensusMEDissimilarity(unmergedMEs);
# Cluster consensus modules
consMETree = hclust(as.dist(consMEDiss), method = "average");
# Plot the result
sizeGrWindow(7,6)
par(mfrow = c(1,1))
plot(consMETree, main = "Consensus clustering of consensus module eigengenes",
     xlab = "", sub = "")
abline(h=0.3, col = "red")
merge = mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.3, verbose = 3)
# Numeric module labels
moduleLabels = merge$colors;
# Convert labels to colors
moduleColors = labels2colors(moduleLabels)
# Eigengenes of the new merged modules:
consMEs = merge$newMEs;

#                         PLOT DENDROGRAM WITH COLOURS
#                   REPRESENTING CO-EXPRESSION CLUSTERS/MODULES
#
#===============================================================================
sizeGrWindow(9,6)
plotDendroAndColors(consTree, cbind(unmergedColors, moduleColors),
                    c("Unmerged", "Merged"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

