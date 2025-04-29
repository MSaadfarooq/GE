# Load libraries ----
suppressPackageStartupMessages({
  library(DESeq2)
  library(vsn)
  library(tidyverse)
  library(tximeta)
  library(apeglm)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(factoextra)
  library(rafalib)
  library(dendextend)
  library(ggdendro)
  library(corrplot)
  library(PCAtools)
})
library(ExploreModelMatrix)
library(hexbin)
library(iSEE)
library(clustertend)
library(NbClust)
library(ClassDiscovery)

# import the sample sheet ----
meta<-read_csv('data/Metadata.csv') # missing IDs ?
#Meta$sampleID <- paste0(meta$cordBlood_ID,'_', meta$bw)

# Set directory path
dir <- "~/results/salmon"
# List directories and extract sample names
dirs <- list.files(dir, full.names = FALSE, pattern = "\\.salmon$")
sample_names <- sub("\\.salmon$", "", dirs)

## tximeta expects a table with at least 2 columns (names and files)
# Modify the file path construction to include the .salmon suffix
coldata <- meta %>%
  select(names = cordBlood_ID, bw, sex, delivery) %>%
  mutate(files = file.path("~/results/salmon",
                           paste0(meta$cordBlood_ID, ".salmon"),"quant.sf"))

# check the modified file paths & if it exists
coldata$files
file.exists(coldata$files)

coldata <- coldata[file.exists(coldata$files),]
stopifnot(all(file.exists(coldata$files)))
# Read in Salmon counts with tximeta to create SE object
se <- tximeta(coldata)
se
# summarise per-transcript counts into per-gene counts
gse <- summarizeToGene(se, assignRanges="abundant")
gse

# number of genes with counts above cutoff 
length(which(rowSums(assay(gse, "counts")) > 10))

# filter for low counts
gse <- gse[rowSums(assay(gse, "counts")) > 10, ]
gse # colData(), rowRanges(), assay()

# assess types of RNA
table(rowData(gse)$gene_biotype)

# subset biotypes — be selective, Avoid very low-confidence types
# keep_types <- c("protein_coding", "lncRNA")
# gse <- gse[rowData(gse)$gene_biotype %in% keep_types, ]

# compare the library sizes of samples
gse$libSize <-  colSums(assay(gse))
colData(gse) |>
  as.data.frame() |>
  ggplot(aes(x = names, y = libSize / 1e6, fill = bw)) + 
  geom_bar(stat = "identity") + theme_bw() + 
  labs(x = "Sample", y = "Total count in millions") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# Write the counts to an object
gene_counts <- assay(gse) %>% 
  round() %>% 
  data.frame()
write.csv(gene_counts, 'output/gene_counts.csv', row.names = TRUE)
save(gse, gene_counts, meta, file = '../cb.rda')  # gse <- load("cb.rda")

# RNA-seq counts distribution
ggplot(gene_counts) +
  geom_histogram(aes(x = CB0003), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

# Raw counts range
summary(gene_counts)

# few outliers affect distribution visualization
boxplot(gene_counts, main='Raw counts', las=2)

# Rows with few counts, if can be ignored
rs <- rowSums(assay(gse))
hist(log10(rs + 1))
abline(v=1, col="blue", lwd=3)

cts <- gene_counts[rs > 100,]

# two samples’ GE values against each other in a scatterplot
smoothScatter(cts[,1:2]) # better job with so many data points

# take the log of counts & plot again
logcts <- log10(cts + 1)
smoothScatter(logcts[,1:2])

avg_expr <- rowMeans(gene_counts)

layout(matrix(1:2, nrow=1))
hist(avg_expr)
hist(log10(avg_expr + 1))

# also log-transform the y-axis
ggplot(data.frame(avg_expr), aes(x=avg_expr)) +
  geom_histogram(bins = 50) +
  scale_x_continuous(breaks = c(0,1,10,100,1000,10000,20000), trans="log1p", expand=c(0,0)) +
  scale_y_continuous(breaks = c(0,1), expand=c(0,0), trans="log1p") +
  theme_bw()

# log2 of counts for transformation
logcounts <- log2(gene_counts + 1)

# make a colour vector
statusCols <- case_when(gse$bw=="low" ~ "red3", 
                        gse$bw=="normal" ~ "orange")
# Check distributions of samples using boxplots
boxplot(logcounts,
        xlab="",
        ylab="Log2(Counts)",
        las=2,
        col=statusCols,
        main="Log2(Counts)")
legend("topright", legend = c("Low", "Normal"),fill = c("red3", "orange"), cex=0.7)

logcounts <- log(gene_counts[,3],10) 
boxplot(logcounts, main="", xlab="", ylab="Raw read counts per gene (log10)",axes=FALSE)
d <- density(logcounts)
plot(d,xlim=c(1,8),main="",ylim=c(0,.45),xlab="Raw read counts per gene (log10)", ylab="Density")

# Subset the read counts for the 100 highly expressed genes ----
top.exp <- order(rowMeans(gene_counts), decreasing=TRUE)[1:100]
top.exp <- gene_counts[top.exp,]

# Annotate the samples in the subset with their age (check order with design!)
colnames(top.exp)<- meta$gestational_age
head(top.exp)
pheatmap(top.exp, cluster_cols = F, labels_row = F)

# compare raw counts per sample, across all genes
sampleCounts <- data.frame(ID=rownames(gse),
                           count=c(assay(gse)),
                           batch=rep(gse$batch, each=nrow(gse)),
                           sex=rep(gse$bw, each=nrow(gse)),
                           condition=rep(gse$bw, each=nrow(gse)))
# draw the ridgeplot
ggplot(sampleCounts %>% filter(count>0),
       aes(x=count, y=condition, fill=sex)) + 
  ggridges::geom_density_ridges(quantile_lines=TRUE, quantiles=2, vline_color="darkred",
                      scale=0.9, alpha=0.5) +
  guides(scale="none") + scale_x_continuous(trans="log10") + theme_bw()

## Create the DESeq dataset object
dds <- DESeqDataSet(gse, design = ~ bw + sex)
dds
# variance increases with the average read count
meanSdPlot(assay(dds), ranks = FALSE)

# visualize the design ----
vd <- VisualizeDesign(sampleData = colData(dds), 
                      designFormula = ~ bw + sex
)
vd
vd$plotlist
# visualize this design
vd <- VisualizeDesign(sampleData = coldata[, c("sex", "bw")], 
                      designMatrix = attr(dds, "modelMatrix"), 
                      flipCoordFitted = TRUE)

# Transform counts for data visualization
vsd <- vst(dds, blind=TRUE)
meanSdPlot(assay(vsd), ranks = FALSE)

# Plot PCA ----
plotPCA(vsd, intgroup="bw") #, ntop = 2000 

# alternate to show variances
pcaData <- plotPCA(vsd, intgroup="sex", returnData=T) 

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = bw), size = 3) +
  theme_bw() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() 

pcaData_libSize <- plotPCA(vsd, intgroup = c("libSize"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData_libSize, "percentVar"))

ggplot(pcaData_libSize, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = libSize / 1e6), size = 3) +
  theme_bw() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()

# Explore the additional PCs or if identify genes that contribute most to the PCs
vsd_mat <- assay(vsd) # Input is a matrix of log transformed values
pcData_add <- prcomp(t(vsd_mat))

# Create data frame with metadata and additional PC values for input to ggplot
df_PC <- cbind(meta, pcData_add$x)
ggplot(df_PC) +
  geom_point(aes(x=PC1, y=PC4, color = bw)) + theme_bw()

# count density for all samples, by group
as.data.frame(vsd_mat) %>%
  pivot_longer(names_to = "sample",
    values_to = "logCounts", everything()) %>% 
  left_join(meta, by = c("sample" = "cordBlood_ID")) %>% 
  ggplot(aes(x=logCounts, group = sample)) +
  geom_density(aes(colour = bw)) +
  scale_colour_manual(values = statusCols) +
  labs(x = "log of Counts", title = "Gene expression count density")

# batches in the PCA plot after vst ----
mm <- model.matrix(~bw, colData(vsd))
mat <- limma::removeBatchEffect(vsd_mat, batch=vsd$batch, design=mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup='bw')

# PCAtools
p <- pca(vsd_mat, metadata = colData(dds), removeVar = 0.1)
screeplot(p, axisLabSize = 12, titleLabSize = 12)
biplot(p, colby = 'bw', shape = 'sex')
biplot(p, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
pairsplot(p)
plotloadings(p, components = getComponents(p, c(1:4)))

# CLUSTERING OF NORMALISED, FILTERED DATA ----
meta$cl <- as.factor(meta$bw)
meta$cl <- factor(meta$cl, levels=c("normal", "low"))
levels(meta$cl) <- c("lightblue","darkblue")
meta$cl <- as.character(meta$cl)

d <- dist(t(vsd_mat))
hc <- hclust(d, method="complete")

# Plot clustering, identifying sample and color coding by status:
pdf("output/NormalisedFiltered_Clustering.pdf", width=25, height=10)
myplclust(hc, labels=meta$cordBlood_ID, lab.col=meta$cl, cex=1.5, main="Samples Clustering (complete)")
legend("topright",legend=c("normal", "low"), cex = 1,
       text.col=c("lightblue","darkblue"), pch=rep(16,2),col=c("lightblue","darkblue"))
dev.off()

# Hierarchical Clustering
hclDat <-  t(vsd_mat) %>% dist(method = "euclidean") %>%
  hclust()

dendro.dat <- as.dendrogram(hclDat) %>% dendro_data()
dendro.dat$labels <- dendro.dat$labels %>%
  left_join(meta, by = c(label = "cordBlood_ID"))

ggplot(dendro.dat$segment) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_label(data = dendro.dat$labels,
             aes(x = x, y = y, label = label,fill = bw),hjust = 0, nudge_y = 1) +
  scale_fill_manual(values = statusCols) +
  coord_flip() +
  labs(x = NULL, y = "Distance", title = NULL, fill = "Sample Group") +
  scale_y_reverse(expand = c(0.3, 0)) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), panel.background = element_blank())

# Cluster Plots 
bw.bar <- c("grey75","darkblue")
bw.bar <- bw.bar[as.numeric(as.factor(meta$bw))]

d <- dist(t(vsd_mat))
hc <- hclust(d, method="average")
dend <- as.dendrogram(hc)
dend <- color_branches(dend, k=2, col=c("darkgreen","darkred"))

# Draw dendrogram
plot(dend)
colored_bars(colors=bw.bar, dend=dend, rowLabels=c("bw")) # needed ?
legend("topright", legend=c("normal", "low"), pch=15, bty="n", col=c("darkblue","grey75"))

# Gender Clustering ----
# X: XIST (ENSG00000229807), Y: RPS4Y1 (ENSG00000129824), Y: EIF1AY (ENSG00000198692)
# Y: DDX3Y (ENSG00000067048), Y: KDM5D (ENSG00000012817)

# Subset the normalised counts to the sex-linked genes:
sex.genes <- c("ENSG00000229807","ENSG00000129824","ENSG00000198692","ENSG00000067048","ENSG00000012817")
norm.sex <- subset(normalized_counts, rownames(normalized_counts) %in% sex.genes)
dim(norm.sex)

# Cluster plot of sex-linked genes:
BW <- bw.bar[as.numeric(as.factor(meta$bw))]
Sex <- c("darkred","orange")
sex_labels <- as.factor(meta$sex)
Sex_colors <- Sex[as.numeric(sex_labels)]

d <- dist(t(norm.sex))
hc <- hclust(d, method="average")
dend <- as.dendrogram(hc)

png("figures/ClustEuclideanAvg_SexGenes.png", width=2800, height=1000, res = 180)
par(mar=c(10,7,2,4))
plot(dend)
colored_bars(colors=Sex_colors, dend=dend, rowLabels="sex")
legend("topright", legend = levels(sex_labels), fill = c("darkred", "orange"), title = "Sex")
dev.off()

# Heatmap + Dendrogram ----
raw.dist <- dist(t(assay(vsd)))
colz <- colorRampPalette(brewer.pal(9, "Blues"))(255)

Heatmap(as.matrix(raw.dist), col = colz,
  name = "Euclidean\ndistance",
  cluster_columns = hclust(raw.dist), show_column_names = FALSE,
  cluster_rows = hclust(raw.dist),
  bottom_annotation = columnAnnotation(Birthweight = vsd$bw, Sex=vsd$sex,
    col = list(Sex = c(female = "pink2", male = "steelblue"),
               Birthweight  = c(low = "forestgreen", normal = "lightblue3")))
) # rearrange similar samples if cluster not required

# Clustering genes most variable across samples ----
VarGenes <- apply(gene_counts, 1, var) # compute variance of each gene

#sort the results by variance in decreasing order and select the top 100 genes
topVarGenes <- names(VarGenes[order(VarGenes, decreasing = T)][1:100])

# overlay some annotation tracks to observer replicate clustering
heatmapColAnnot <- data.frame(colData(gse)[, c("bw", "sex")])
pheatmap(gene_counts[topVarGenes,], scale = 'row',
         show_rownames = FALSE,
         annotation_col = heatmapColAnnot)

# Correlation plot ----
CorMat <- cor(vsd_mat) # create cor matrix

corrplot(CorMat, order = 'hclust',
         addrect = 2, addCoef.col = 'white',  # split clusters into group & surround with rectangle
         number.cex = 0.7)
corrplot(CorMat, method = 'square', type = 'lower', diag = FALSE)

# split the clusters into *2* based on the clustering similarity
pheatmap(CorMat, # border_color=NA,
         color = colz, annotation_col = heatmapColAnnot,
         cutree_cols = 2) # may use cutree_rows

# Corplot of top var genes in all samples
corVar <- cor(vsd_mat[topVarGenes,])
rownames(corVar) <- str_c(meta$bw, "_", meta$sex)
colnames(corVar) <- str_c(meta$bw, "_", meta$sex)
corrplot(corVar, 
         method = "color", 
         addCoef.col = "black", 
         number.digits= 3,
         order = "hclust",
         is.corr = FALSE)

# Highly expressed genes common in PL & CB groups ----
PL_geneTop <- as.data.frame(top.exp) %>%
  rownames_to_column("Gene")
CB_geneTop <- as.data.frame(top.expCB) %>%
  rownames_to_column("Gene")

genesTop_PL_CB <- inner_join(PL_geneTop, CB_geneTop, by = "Gene")

# add symbols to them after next session
genesTop_PL_CB <- genesTop_PL_CB %>%
  inner_join(anno.genes, by = c("Gene" = "gene_id"))
write.csv(genesTop_PL_CB, "Hi_Exp_CB-PL.csv")

TopExpVenn <- list(Placenta=PL_geneTop$Gene, CordBlood=CB_geneTop$Gene)

# draw Venn Diagram ----
ggvenn(TopExpVenn, fill_color = c("lightblue3", "magenta4"), stroke_color = F)

# Calculate mean expression across samples for each gene ----
mean_pl <- rowMeans(PL_gene_counts)
mean_cb <- rowMeans(gene_counts)

# Combine the data into a single data frame
df_comparison <- tibble(
  Mean_Expression = c(mean_pl, mean_cb),
  Tissue = c(rep("pl", length(mean_pl)), rep("cb", length(mean_cb)))
)
# Apply log2 transformation to the data (adding 1 to avoid log(0))
df_comparison <- df_comparison %>%
  mutate(Mean_Expression_Log = log2(Mean_Expression + 1))

# Violin plot to visualize the distribution of the data
ggplot(df_comparison, aes(x = Condition, y = Mean_Expression_Log, fill = Condition)) +
  geom_violin(trim = FALSE) +  # Trim = FALSE includes the full distribution
  labs(title = "Log-Transformed Gene Expression Comparison between pl and cb",
       y = "Mean Expression", x = "Tissue") +
  scale_fill_manual(values = c("pl" = "royalblue", "cb" = "darkred"))

# Bar plot with mean and error bars (using standard deviation)
ggplot(df_comparison, aes(x = Tissue, y = Mean_Expression_Log, fill = Tissue)) +
  stat_summary(fun = "mean", geom = "bar", position = "dodge", width = 0.6) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2) +
  labs(title = "Log-Transformed Gene Expression Comparison between pl and cb",
       y = "Mean Expression", x = "Tissue") +
  scale_fill_manual(values = c("pl" = "royalblue", "cb" = "darkred"))

# Density plot to compare distributions between tissues
ggplot(df_comparison, aes(x = Mean_Expression_Log, fill = Tissue)) +
  geom_density(alpha = 0.4) +  # Use alpha to adjust transparency
  labs(title = "Log-Transformed Gene Expression Comparison between pl and cb",
       y = "Mean Expression", x = "Tissue") + theme_bw()+
  scale_fill_manual(values = c("pl" = "royalblue", "cb" = "darkred"))
