# Load libraries ----
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(vsn)
  library(tidyverse)
  library(tximeta)
  library(apeglm)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(hexbin)
  library(iSEE)
  library(ClassDiscovery)
  library(factoextra)
  library(clustertend)
  library(NbClust)
  library(rafalib)
  library(dendextend)
  library(ggdendro)
  library(PCAtools)
})

# import the sample sheet ----
meta<-read_csv('Metadata.csv')

# subset meta for placenta, Set directory path
dir <- "/home/saad/rna_seq/EDA/salmon/"
# List directories and extract sample names
dirs <- list.files(dir, full.names = FALSE, pattern = "\\.salmon$")
sample_names <- sub("\\.salmon$", "", dirs)

# Filter the dataframe
Meta_PL <- meta[meta$placenta_ID %in% sample_names, ]

## tximeta expects a table with at least 2 columns (names and files)
# Modify the file path construction to include the .salmon suffix
coldata <- Meta_PL %>%
  dplyr::select(names = placenta_ID, bw, sex) %>%
  mutate(files = file.path("salmon", paste0(Meta_PL$placenta_ID, ".salmon"), "quant.sf"),
         bw = factor(bw), bw = relevel(bw, 'low'))

# Check the modified file paths
coldata$files

# Check if the files exist now
file.exists(coldata$files)

# Read in Salmon counts with tximeta to create SE object
se <- tximeta(coldata)
se
# summarise per-transcript counts into per-gene counts
PL_gse<-summarizeToGene(se, assignRanges="abundant")
PL_gse
# filter for low counts
PL_gse <- PL_gse[rowSums(assay(PL_gse, "counts")) > 5, ]
PL_gse # colData(), rowRanges(), assay()

# number of genes with counts above 50 i.e., TPM > 1
length(which(rowSums(assay(PL_gse, "counts")) > 50))

# assess types of RNA
table(rowData(PL_gse)$gene_biotype)

# subsetting to only the mRNA genes
pl_mRNA <- PL_gse[rowData(PL_gse)$gene_biotype == "protein_coding", ]
dim(pl_mRNA)
pl_mRNA@rowRanges@ranges@NAMES |> head()

length(which(rowSums(assay(pl_mRNA, "counts")) > 50))

# compare the library sizes of samples
PL_gse$libSize <-  colSums(assay(PL_gse))
colData(PL_gse) |>
  as.data.frame() |>
  ggplot(aes(x = names, y = libSize / 1e6, fill = bw)) + 
  geom_bar(stat = "identity") + theme_bw() + 
  labs(x = "Sample", y = "Total count in millions") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# Write the counts to an object
PL_gene_counts <- assay(PL_gse) %>% 
  round() %>% 
  data.frame()
saveRDS(PL_gse, 'PL_gse.RDS')
write.csv(PL_gene_counts, 'PL_gene_counts.csv', row.names = TRUE)

# RNA-seq counts distribution
ggplot(PL_gene_counts) +
  geom_histogram(aes(x = PL0007), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

# Raw counts range
summary(PL_gene_counts)

# few outliers affect distribution visualization
boxplot(PL_gene_counts, main='Raw counts', las=2)

# Rows with few counts, if can be ignored
rs <- rowSums(assay(PL_gse))
hist(log10(rs + 1))
abline(v=1, col="blue", lwd=3)

cts <- PL_gene_counts[rs > 100,]

# two samples’ GE values against each other in a scatterplot
smoothScatter(cts[,1:2]) # better job with so many data points

# take the log of counts & plot again
logcts <- log10(cts + 1)
smoothScatter(logcts[,1:2])

# get a sense for the distribution of counts (log-scale)
hist(logcts[,1])

avg_expr <- rowMeans(PL_gene_counts)

layout(matrix(1:2, nrow=1))
hist(avg_expr)
hist(log10(avg_expr + 1))

# also log-transform the y-axis
ggplot(data.frame(avg_expr), aes(x=avg_expr)) +
  geom_histogram(bins = 50) +
  scale_x_continuous(breaks = c(0,1,10,100,1000,10000,20000), trans="log1p", expand=c(0,0)) +
  scale_y_continuous(breaks = c(0,1), expand=c(0,0), trans="log1p")

# log2 of counts for transformation
logcounts <- log2(PL_gene_counts + 1)

# make a colour vector
statusCols <- case_when(PL_gse$bw=="low" ~ "red3", 
                        PL_gse$bw=="normal" ~ "orange")
# Check distributions of samples using boxplots
boxplot(logcounts,
        xlab="",
        ylab="Log2(Counts)",
        las=2,
        col=statusCols,
        main="Log2(Counts)")
legend("topright", legend = c("Low", "Normal"),fill = c("red3", "orange"), cex=0.7)

logcounts <- log(PL_gene_counts[,3],10) 
boxplot(logcounts, main="", xlab="", ylab="Raw read counts per gene (log10)",axes=FALSE)
d <- density(logcounts)
plot(d,xlim=c(1,8),main="",ylim=c(0,.45),xlab="Raw read counts per gene (log10)", ylab="Density")

# Subset the read counts object for the 100 most highly expressed genes
top.exp <- order(rowMeans(PL_gene_counts), decreasing=TRUE)[1:100]
top.exp <- PL_gene_counts[top.exp,]

# Annotate the samples in the subset with their age (check order with design!)
colnames(top.exp)<- Meta_PL$gestational_age
head(top.exp)
pheatmap(top.exp, cluster_cols = F, labels_row = F)

## Create the DESeq dataset object
dds <- DESeqDataSet(PL_gse, design = ~ sex + bw)
dds
levels(dds$bw) # dds$Condition <- relevel(dds$Condition, ref = "mock") to relevel

# variance increases with the average read count
meanSdPlot(assay(dds), ranks = FALSE)

# visualize this design
vd <- VisualizeDesign(sampleData = colData(dds)[, c("sex", "bw")], 
                      designMatrix = attr(dds, "modelMatrix"), 
                      flipCoordFitted = TRUE)
vd

# Transform counts for data visualization ----
PL_rld <- vst(dds, blind=TRUE)
meanSdPlot(assay(PL_rld), ranks = FALSE)

# Plot PCA, ntop= can be used 
plotPCA(PL_rld, intgroup="bw")

# alternate to show variances
pcaData <- plotPCA(PL_rld, intgroup="sex", returnData=T)

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = sex), size = 3) +
  theme_bw() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() #scale_color_manual(values = c(Male = "royalblue", Female = "red3"))

pcaData_PL_rld <- DESeq2::plotPCA(PL_rld, intgroup = c("libSize"),
                              returnData = TRUE)
percentVar <- round(100 * attr(pcaData_PL_rld, "percentVar"))

ggplot(pcaData_PL_rld, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = libSize / 1e6), size = 3) +
  theme_bw() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()

# Explore the additional PCs in your data or if you would like to identify genes that contribute most to the PCs. Input is a matrix of log transformed values
PL_rld_mat <- assay(PL_rld)
pca <- prcomp(t(PL_rld_mat))

# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(Meta_PL, pca$x)
ggplot(df) +
  geom_point(aes(x=PC1, y=PC3, color = bw)) + theme_bw()

# CLUSTERING OF NORMALISED, FILTERED DATA ----
Meta_PL$cl_bw <- as.factor(Meta_PL$bw)
Meta_PL$cl_bw <- factor(Meta_PL$cl_bw, levels=c("normal", "low"))
levels(Meta_PL$cl_bw) <- c("lightblue","darkblue")
Meta_PL$cl_bw <- as.character(Meta_PL$cl_bw)

d <- dist(t(PL_rld_mat))
hc <- hclust(d, method="complete")

# Plot clustering, identifying sample and color coding by status:
pdf("PL/NormalisedFiltered_Clustering_bw.pdf", width=25, height=10)
myplclust(hc, labels=Meta_PL$placenta_ID, lab.col=Meta_PL$cl_bw, cex=1.5, main="Placenta Samples Clustering (complete)")
legend("topright",legend=c("Normal","Low"), cex = 1,
       text.col=c("lightblue","darkblue"), pch=rep(16,2),col=c("lightblue","darkblue"))
dev.off()

# Hierarchical Clustering ----
hclDat <-  t(PL_rld_mat) %>%
  dist(method = "euclidean") %>%
  hclust()
ggdendrogram(hclDat, rotate=TRUE, scale.color=statusCols)

# CLUSTER PLOTS 
bw.bar <- c("grey75","darkblue")
bw.bar <- bw.bar[as.numeric(as.factor(Meta_PL$bw))]

d <- dist(t(PL_rld_mat))
hc <- hclust(d, method="average")
dend <- as.dendrogram(hc)
dend <- color_branches(dend, k=2, col=c("darkgreen","darkred"))

# Draw dendrogram
plot(dend)
colored_bars(colors=bw.bar, dend=dend, rowLabels=c("bw")) # needed ?
legend("topright", legend=c("normal","low"), pch=15, bty="n", col=c("darkblue","grey75"))

# CHECK GENDER CLUSTERING ----
################

# X: XIST (ENSG00000229807)
# Y: RPS4Y1 (ENSG00000129824)
# Y: EIF1AY (ENSG00000198692)
# Y: DDX3Y (ENSG00000067048)
# Y: KDM5D (ENSG00000012817)

# Subset the normalised counts to the sex-linked genes above:
sex.genes <- c("ENSG00000229807","ENSG00000129824","ENSG00000198692","ENSG00000067048","ENSG00000012817")
norm.sex <- subset(PL_normalized_counts, rownames(PL_normalized_counts) %in% sex.genes)

dim(norm.sex)

# Cluster plot of sex-linked genes:
BW <- c("skyblue","grey75")
BW <- BW[as.numeric(as.factor(Meta_PL$bw))]
Sex <- c("darkred","darkblue")
Sex <- Sex[as.numeric(as.factor(Meta_PL$newborn_gender))]

d <- dist(t(norm.sex))
hc <- hclust(d, method="average")
dend <- as.dendrogram(hc)

pdf("EDA_PL/SexLinkGenes_EuclideanAvgClustering.pdf", width=50, height=10)
par(mar=c(10,7,2,4))
plot(dend)
colored_bars(colors=cbind(BW,Sex), dend=dend, rowLabels=c("BW","Sex"))
dev.off()

# Dendrogram + Heatmap to visualize Euclidean distances ----
raw.dist <- dist(t(assay(PL_rld)))

colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)

ComplexHeatmap::Heatmap(
  as.matrix(raw.dist), 
  col = colors,
  name = "Euclidean\ndistance",
  cluster_rows = hclust(raw.dist),
  cluster_columns = hclust(raw.dist),
  bottom_annotation = columnAnnotation(
    birthweight = PL_rld$bw ,
    sex = PL_rld$sex,
    col = list(sex = c(Female = "pink2", Male = "lightblue3"),
               birthweight  = c(low = "forestgreen", normal = "purple3")))
) 

#https://tavareshugo.github.io/data-carpentry-rnaseq/02_rnaseq_exploratory.html
cts <- PL_gene_counts %>% 
  pivot_longer(cols = PL0002:PL0080, 
               names_to = "sample", 
               values_to = "cts")
cts %>% 
  pivot_wider(names_from = "sample", values_from = "cts")

Meta_PL$sample <- Meta_PL$placenta_ID
trans_cts_long <- full_join(cts, Meta_PL, by = 'sample')

trans_cts_long %>%
  ggplot(aes(cts, colour = sex)) + 
  geom_freqpoly() + 
  facet_grid(bw ~ gestational_age)

# Covariation between between samples
PL_gene_counts %>%  # wide format of counts table
  ggplot(aes(PL0002, PL0009)) + geom_point() +
  geom_abline(colour = "brown")