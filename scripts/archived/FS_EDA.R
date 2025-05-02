# Load libraries ----
suppressPackageStartupMessages({
  library(DESeq2)
  library(vsn)
  library(tidyverse)
  library(tximeta)
  library(apeglm)
  library(ComplexHeatmap)
  library(pheatmap)
  library(RColorBrewer)
  library(factoextra)
  library(clustertend)
  library(NbClust)
  library(rafalib)
  library(dendextend)
  library(ggdendro)
  library(corrplot)
})

# import the sample sheet ----
meta<-read_csv('Metadata.csv')

# subset meta for placenta, Set directory path
dir <- "/home/saad/rna_seq/fastq_done/results_salmon/"
# List directories and extract sample names
dirs <- list.files(dir, full.names = FALSE, pattern = "\\.salmon$")
sample_names <- sub("\\.salmon$", "", dirs)

# Filter the dataframe
meta_fs <- meta[meta$fatherSemen_ID %in% sample_names, ]

## tximeta expects a table with at least 2 columns (names and files)
# Modify the file path construction to include the .salmon suffix
coldata <- meta_fs %>%
  dplyr::select(names = fatherSemen_ID, bw, sex) %>%
  mutate(files = file.path("../fastq_done/results_salmon",
                           paste0(meta_fs$fatherSemen_ID, ".salmon"),"quant.sf"))

# Check the modified file paths
coldata$files

# Check if the files exist now
file.exists(coldata$files)

# Read in Salmon counts with tximeta to create SE object
se <- tximeta(coldata)
se
# summarise per-transcript counts into per-gene counts
gse<-summarizeToGene(se, assignRanges="abundant")
gse
# filter for low counts
gse <- gse[rowSums(assay(gse, "counts")) > 10, ]
gse # colData(), rowRanges(), assay()

# number of genes with counts above 50 i.e., TPM > 1
length(which(rowSums(assay(gse, "counts")) > 50))

# assess types of RNA
table(rowData(gse)$gene_biotype)

# subsetting to only the mRNA genes


# compare the library sizes of samples
gse$libSize <-  colSums(assay(gse))
colData(gse) |>
  as.data.frame() |>
  ggplot(aes(x = names, y = libSize / 1e6, fill = bw)) + 
  geom_bar(stat = "identity") + theme_bw() + 
  labs(x = "Sample", y = "Total count in millions") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave('FS/EDA/libsize.png')
# Write the counts to an object
gene_counts <- assay(gse) %>% 
  round() %>% 
  data.frame()
saveRDS(gse, 'FS/EDA/gse.RDS')
write.csv(gene_counts, 'FS/EDA/gene_counts.csv', row.names = TRUE)

# RNA-seq counts distribution
#ggplot(gene_counts) +
 # geom_histogram(aes(x = CB0007), stat = "bin", bins = 200) +
  #xlab("Raw expression counts") +
  #ylab("Number of genes")

# Raw counts range
#summary(gene_counts)

# few outliers affect distribution visualization
#boxplot(gene_counts, main='Raw counts', las=2)


# log2 of counts for transformation
logcounts <- log2(gene_counts + 1)

# make a colour vector
statusCols <- case_when(gse$bw=="low" ~ "red3", 
                        gse$bw=="normal" ~ "orange")
# Check distributions of samples using boxplots
pdf('FS/EDA/count_dist.pdf')
boxplot(logcounts,
        xlab="",
        ylab="Log2(Counts)",
        las=2,
        col=statusCols,
        main="Log2(Counts)")
legend("topright", legend = c("Low", "Normal"),fill = c("red3", "orange"), cex=0.7)
dev.off()

# Subset the read counts for the 100 highly expressed genes ----
top.exp <- order(rowMeans(gene_counts), decreasing=TRUE)[1:100]
top.exp <- gene_counts[top.exp,]

# Annotate the samples in the subset with their age (check order with design!)
colnames(top.exp)<- meta_fs$gestational_age
head(top.exp)
pdf('FS/EDA/hiExp100.pdf')
pheatmap(top.exp, cluster_cols = F, labels_row = F)
dev.off()

# Clustering genes most variable across samples ----
VarGenes <- apply(gene_counts, 1, var) # compute variance of each gene

#sort the results by variance in decreasing order and select the top 100 genes
topVarGenes <- names(VarGenes[order(VarGenes, decreasing = T)][1:100])
pdf('FS/EDA/hiVar100.pdf')
pheatmap(gene_counts[topVarGenes,], scale = 'row', show_rownames = FALSE)
dev.off()
# overlay some annotation tracks to observer replicate clustering
heatmapColAnnot <- data.frame(colData(gse)[, c("bw", "sex")])
#pheatmap(gene_counts[topVarGenes,], scale = 'row',
         #show_rownames = FALSE,
         #annotation_col = heatmapColAnnot)

## Create the DESeq dataset object
dds <- DESeqDataSet(gse, design = ~ bw)
dds
# variance increases with the average read count
#meanSdPlot(assay(dds), ranks = FALSE)

# visualize this design
#vd <- VisualizeDesign(sampleData = colData(dds)[, c("sex", "bw")], 
                      #designMatrix = attr(dds, "modelMatrix"), 
                      #flipCoordFitted = TRUE)
#vd

# Transform counts for data visualization ----
vsd <- vst(dds, blind=TRUE)
#meanSdPlot(assay(vsd), ranks = FALSE)

# Plot PCA, ntop= can be used
pdf('FS/EDA/pca_bw.pdf')
plotPCA(vsd, intgroup="bw")
dev.off()
# alternate to show variances
pcaData <- plotPCA(vsd, intgroup="sex", returnData=T)

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = sex), size = 3) +
  theme_bw() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() 
ggsave('FS/EDA/pca_bw_sex.png')
pcaData_vsd <- DESeq2::plotPCA(vsd, intgroup = c("libSize"),
                              returnData = TRUE)
percentVar <- round(100 * attr(pcaData_vsd, "percentVar"))

ggplot(pcaData_vsd, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = libSize / 1e6), size = 3) +
  theme_bw() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
ggsave('FS/EDA/pca_lib.png')
# Explore the additional PCs in your data or if you would like to identify genes that contribute most to the PCs. Input is a matrix of log transformed values
vsd_mat <- assay(vsd)
pcData <- prcomp(t(vsd_mat))

# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(meta_fs, pcData$x)
ggplot(df) +
  geom_point(aes(x=PC1, y=PC3, color = bw)) + theme_bw()
ggsave('FS/EDA/pc3-4.png')
# PCA plot with sample lables 
#autoplot(pcData, data = meta_fs, colour = 'bw', shape='sex', size=3) +
 # geom_text_repel(aes(x = PC1, y = PC2, label = cordBlood_ID), box.padding = 0.8)

# CLUSTERING OF NORMALISED, FILTERED DATA ----
meta_fs$cl_bw <- as.factor(meta_fs$bw)
meta_fs$cl_bw <- factor(meta_fs$cl_bw, levels=c("normal", "low"))
levels(meta_fs$cl_bw) <- c("lightblue","darkblue")
meta_fs$cl_bw <- as.character(meta_fs$cl_bw)

d <- dist(t(vsd_mat))
hc <- hclust(d, method="complete")

# Plot clustering, identifying sample and color coding by status:
pdf("FS/EDA/NormalisedFiltered_Cluster_bw.pdf", width=25, height=10)
myplclust(hc, labels=meta_fs$cordBlood_ID, lab.col=meta_fs$cl_bw, cex=1.5, main="CordBlood Samples Clustering (complete)")
legend("topright",legend=c("Normal","Low"), cex = 1,
       text.col=c("lightblue","darkblue"), pch=rep(16,2),col=c("lightblue","darkblue"))
dev.off()

# Hierarchical Clustering ----
hclDat <-  t(vsd_mat) %>%
  dist(method = "euclidean") %>%
  hclust()
#ggdendrogram(hclDat, rotate=TRUE, scale.color=statusCols)

# more customized dendrogram


# CLUSTER PLOTS 
bw.bar <- c("grey75","darkblue")
bw.bar <- bw.bar[as.numeric(as.factor(meta_fs$bw))]

d <- dist(t(vsd_mat))
hc <- hclust(d, method="average")
dend <- as.dendrogram(hc)
dend <- color_branches(dend, k=2, col=c("darkgreen","darkred"))

# Draw dendrogram
pdf('FS/EDA/dendro.pdf')
plot(dend)
colored_bars(colors=bw.bar, dend=dend, rowLabels=c("bw")) # needed ?
legend("topright", legend=c("normal","low"), pch=15, bty="n", col=c("darkblue","grey75"))
dev.off()

# Dendrogram + Heatmap to visualize Euclidean distances ----
raw.dist <- dist(t(assay(vsd)))

colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)
pdf('FS/EDA/hm_bw_sex.pdf')
ComplexHeatmap::Heatmap(
  as.matrix(raw.dist), 
  col = colors,
  name = "Euclidean\ndistance",
  cluster_rows = hclust(raw.dist),
  cluster_columns = hclust(raw.dist),
  bottom_annotation = columnAnnotation(
    birthweight = vsd$bw ,
    sex = vsd$sex,
    col = list(sex = c(Female = "pink2", Male = "lightblue3"),
               birthweight  = c(low = "forestgreen", normal = "purple3")))
)
dev.off()
# rearrange similar samples if clustering not required

# Correlation plot ----
CorMat <- cor(vsd_mat) # create cor matrix

corrplot(CorMat, order = 'hclust',
         addrect = 2, addCoef.col = 'white',  # split clusters into group & surround with rectangle
         number.cex = 0.7)
corrplot(CorMat, method = 'square', type = 'lower', diag = FALSE)

# clustering + heatmap of correlation matrix
heat.colors <- brewer.pal(6, "Blues")
pdf('FS/EDA/cor_hm.pdf')
pheatmap(CorMat, annotation_col = heatmapColAnnot, color = heat.colors,
         border_color=NA, fontsize = 10, fontsize_row = 10, height=20)
dev.off()


