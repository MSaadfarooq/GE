# Attach the library
library(tidyverse)
library(WGCNA)
library(DESeq2)
library(igraph)
enableWGCNAThreads(nThreads = 24)

# Read in metadata TSV file
Meta_CB <- readr::read_csv('CB/Meta_CB.csv')

# Read in normalized count file after DESeq output
input_mat <- read.csv('FS/DE/normalized_counts.csv')
# input_mat <- column_to_rownames(input_mat, var = "X")

# NOTE: input can be from DESeq -> vst -> normalize counts 
# what if want to work for only 1 group? without model? log2 transformation of the TPM values?
## Alternatively create DESeq dataset object without specifying model after loading 'gse'
#dds <- DESeqDataSet(CB_gse, design = ~ 1)

# Normalize and transform the data in the `DESeqDataSet` object using the `vst()`
#dds <- vst(dds)
# Retrieve the normalized data from the `DESeqDataSet`
#normalized_counts <- assay(dds) 
# normalized_counts <- counts(dds, normalized = T) * not worked

expr_normalized_df <- data.frame(normalized_counts) %>%
  mutate(
    Gene_id = row.names(normalized_counts)
  ) %>%
  pivot_longer(-Gene_id)

expr_normalized_df %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "Normalized Expression",
    x = "Samples",
    y = "Expression values"
  )
input_mat <- t(input_mat)

input_mat[1:5,1:10]           # Look at first 5 rows and 10 columns

# QC of samples with too many missing data as well as genes with zero variance
gsg <- goodSamplesGenes(datExpr = input_mat, verbose = 3)
print(paste("All ok?", gsg$allOK))

# may need to remove problematic samples and feature
sprintf("Removing %d features", ncol(input_mat) - sum(gsg$goodGenes))
sprintf("Removing %d samples", nrow(input_mat) - sum(gsg$goodSamples))
input_mat <- input_mat[gsg$goodSamples, gsg$goodGenes]

# identify outlier samples by hierarchical clustering
sampleTree <- hclust(dist(input_mat), method = "average") #Clustering samples based on distance 

#Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
#draw on line to show cutoff height
abline(h = 95, col = "red"); # h varies depending on plot

# cut height of 15 would cut xx and retain the rest of the samples
cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 15, minSize = 10) # returns numeric vector
input_mat <- input_mat[cut.sampleTree==1, ] # Remove outlier

# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

sft <- pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)
# Power to be used for dataset
sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  geom_point() +
  # put the Power labels slightly above the data points
  geom_text(nudge_y = 0.02) +
  # plot R^2 cutoff
  geom_hline(yintercept = 0.9, col = "red") +
  # Just in case our values are low, we want to make sure we can still see the 0.8 level
  ylim(c(min(sft_df$model_fit), 1.05)) +
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence")
ggsave(filename = 'CB/sft.png')
pdf("CB/sftP.pdf", width=40, height=7)
cex1 = 1
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")
dev.off()
sprintf("Optimal soft-power = %d", sft$powerEstimate)

# visualize the difference between the raw correlations & the soft-thresholded correlations
corRaw <- abs(cor(input_mat))
corSoft <- abs(corRaw**sft$powerEstimate)
par(mfrow=c(1,2))
hist(corRaw, main="Raw correlations")
hist(corSoft, main="Soft-thresholded correlations")

# Pick a soft threshold power near the curve of the Scale independence plot
picked_power <- 9
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
bwnet <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 6000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)
cor <- temp_cor     # Return cor function to original namespace
# save WGCNA results
write_rds(bwnet, file = "Outputs/CB_wgcna_result.RDS")

# Explore WGCNA results ----
# Convert labels to colors for plotting
mergedColors <- labels2colors(bwnet$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  bwnet$dendrograms[[1]],
  mergedColors[bwnet$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05)

bwnet$colors[bwnet$blockGenes[[1]]]
table(bwnet$colors)

# dataframe of eigengene module for each sample in the MEs slot
module_eigengenes <- bwnet$MEs
head(module_eigengenes)

# Module trait assocaition- pull out the list of modules
module_df <- data.frame(
  gene_id = names(bwnet$colors),
  colors = labels2colors(bwnet$colors)
)
module_df[1:5,]
write_delim(module_df,
            file = "FS/gene_modules.txt",
            delim = "\t")
# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order <- names(MEs0) %>% gsub("ME","", .)

# Add treatment names ** sample or Group ???
MEs0$Sample <- row.names(MEs0)

# tidy & plot data
mME <- MEs0 %>%
  pivot_longer(-Sample) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=Sample, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-Trait Relationships", y = "Modules", fill="corr")
ggsave('CB/module-trait-relation.png')
# heatmap of adjacencies among eigengenes & dendrogram of their relationship
plotEigengeneNetworks(MEs0, '', marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))

# Examine Expression Profiles ----
# pick out a few modules of interest here
modules_of_interest <- c("brown1", "cornsilk", "yellowgreen")

# Pull out list of genes in that module
submod <- module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes - use read.csv
normalized_counts[1:5,1:10]

subexpr <- normalized_counts[submod$gene_id,]

submod_df <- data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module), alpha = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(rows = vars(module)) +
  labs(x = "Samples", y = "Normalized Expression")

# Export network file for Cytoscape or as an edge/vertices file
genes_of_interest <- module_df %>%
  subset(colors %in% modules_of_interest)

expr_of_interest <- normalized_counts[genes_of_interest$gene_id,]
expr_of_interest[1:5,1:5]

# Only recalculate TOM for modules of interest
TOM <- TOMsimilarityFromExpr(t(expr_of_interest),
                            power = picked_power)
# Add gene names to row and columns
row.names(TOM) <- row.names(expr_of_interest)
colnames(TOM) <- row.names(expr_of_interest)

edge_list <- data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )

head(edge_list)

# Export Network file to be read into Cytoscape, VisANT, etc
write_delim(edge_list,
            file = "CB/edgelist.tsv",
            delim = "\t")

# Gene relationship to traits and modules ----
metadataNum <- data.frame("Birthweight"=scale(Meta_CB$newborn_weight_kg), row.names = Meta_CB$cordBlood_ID)
nSamples <- nrow(input_mat)
moduleMembership <- as.data.frame(cor(input_mat, MEs0))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(moduleMembership), nSamples))

modNames <- substring(names(MEs0), 3)
names(moduleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

head(moduleMembership[,1:5])
head(MMPvalue[,1:5])

geneTraitSignificance <- as.data.frame(cor(input_mat, metadataNum, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS", names(metadataNum), sep = "")
names(GSPvalue) <- paste("p.GS", names(metadataNum), sep = "")

names(mergedColors) <- colnames(input_mat)
module <- 'pink'
moduleGenes <- names(mergedColors)[which(mergedColors==module)]
geneTraitSignifColumn <- "GSSubtype"
# Make sure that the geneTraitSignificance column is numeric
geneTraitSignificance[, "GSBirthweight"] <- as.numeric(as.character(geneTraitSignificance[, "GSBirthweight"]))

# Generate the scatter plot
verboseScatterplot(
  abs(moduleMembership[moduleGenes, paste('MM', module, sep="")]),
  abs(geneTraitSignificance[moduleGenes, "GSBirthweight"]),
  xlab = paste("Module Membership (MM) in", module, "module"),
  ylab = paste("Gene Significance (GS) for Birthweight"),
  main = paste("Module Membership (MM) vs Gene Significance (GS)\n"),
  cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module
)
resultsModuleMembership <- cbind(moduleMembership, MMPvalue)
mmNames <- names(resultsModuleMembership)
mmNames <- mmNames[order(gsub("p\\.", "", mmNames))]
resultsModuleMembership <- resultsModuleMembership[,mmNames]

resultsSignificance <- cbind(geneTraitSignificance, GSPvalue)
gsNames <- names(resultsSignificance)
gsNames <- gsNames[order(gsub("p\\.", "", gsNames))]
resultsSignificance <- resultsSignificance[,gsNames]


results <- data.frame("Gene"=names(mergedColors), "Module"=unlist(mergedColors), 
                     resultsModuleMembership, resultsSignificance)
results <- results[order(results$Module),]
write.table(results, file="Outputs/Module_results.tsv", sep="\t", row.names=F, col.names=T, quote=F)
head(results[,1:6])

# Modules with biggest differences across groups [refinebio]
# What genes are a part of module 19/or color?

# subset down Network by weight or minimal spanning to identify hub genes
el <- as.data.frame(edge_list)

threshold <- 0.5
filtered_edges <- el[el$correlation > threshold, ]

# Create igraph object
g <- graph_from_data_frame(d = filtered_edges, directed = FALSE)

# Plot the network
plot(g, vertex.label = NA, edge.arrow.size = 0.5)

plot(g, 
     vertex.size = 5, 
     vertex.color = "skyblue", 
     vertex.frame.color = "white", 
     edge.color = "grey", 
     edge.width = edge_list$correlation * 10,  # Scale edge width by correlation
     layout = layout_with_fr)  # Use Fruchterman-Reingold layout

# Which modules have biggest differences across  groups ----
module_eigengenes <- read_csv("CB/module_eigengenes.csv")
module_eigengenes <- column_to_rownames(module_eigengenes, var = '...1')

all.equal(Meta_CB$cordBlood_ID, rownames(module_eigengenes))

# Create the design matrix from the `bw` variable
des_mat <- model.matrix(~ Meta_CB$bw)

# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(module_eigengenes), design = des_mat)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")
head(stats_df)

# plot module that has highest difference
MEmediumspringgreen_df <- module_eigengenes %>%
  rownames_to_column("sampleID") %>%
  # Here we are performing an inner join with a subset of metadata
  inner_join(Meta_CB %>%
               select(cordBlood_ID, bw),
             by = c("sampleID" = "cordBlood_ID")
  )
ggplot(
  MEmediumspringgreen_df,
  aes(
    x = bw,
    y = MEmediumspringgreen,
    color = bw
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()
