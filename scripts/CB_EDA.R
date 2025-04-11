# Load libraries ----
suppressPackageStartupMessages({
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
  library(corrplot)
  library(PCAtools)
})

# import the sample sheet ----
meta<-read_csv('Metadata.csv')

# subset meta for placenta, Set directory path
dir <- "/home/saad/rna_seq/fastq_done/results_salmon/"
# List directories and extract sample names
dirs <- list.files(dir, full.names = FALSE, pattern = "\\.salmon$")
sample_names <- sub("\\.salmon$", "", dirs)

# Filter the dataframe
Meta_CB <- meta[meta$cordBlood_ID %in% sample_names, ]

## tximeta expects a table with at least 2 columns (names and files)
# Modify the file path construction to include the .salmon suffix
coldata <- Meta_CB %>%
  dplyr::select(names = cordBlood_ID, bw, sex) %>%
  mutate(files = file.path("../fastq_done/results_salmon",
                           paste0(Meta_CB$cordBlood_ID, ".salmon"),"quant.sf"))

# Check the modified file paths
coldata$files

# Check if the files exist now
file.exists(coldata$files)

# Read in Salmon counts with tximeta to create SE object
se <- tximeta(coldata)
se
# summarise per-transcript counts into per-gene counts
CB_gse<-summarizeToGene(se, assignRanges="abundant")
CB_gse
# filter for low counts
CB_gse <- CB_gse[rowSums(assay(CB_gse, "counts")) > 10, ]
CB_gse # colData(), rowRanges(), assay()

# number of genes with counts above 50 i.e., TPM > 1
length(which(rowSums(assay(CB_gse, "counts")) > 50))

# assess types of RNA
table(rowData(CB_gse)$gene_biotype)

# subsetting to only the mRNA genes
CB_mRNA <- CB_gse[rowData(CB_gse)$gene_biotype == "protein_coding", ]
dim(CB_mRNA)
CB_mRNA@rowRanges@ranges@NAMES |> head()

length(which(rowSums(assay(CB_mRNA, "counts")) > 20))

# compare the library sizes of samples
CB_gse$libSize <-  colSums(assay(CB_gse))
colData(CB_gse) |>
  as.data.frame() |>
  ggplot(aes(x = names, y = libSize / 1e6, fill = bw)) + 
  geom_bar(stat = "identity") + theme_bw() + 
  labs(x = "Sample", y = "Total count in millions") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# Write the counts to an object
CB_gene_counts <- assay(CB_gse) %>% 
  round() %>% 
  data.frame()
saveRDS(CB_gse, 'EDA_CB/CB_gse.RDS')
write.csv(CB_gene_counts, 'EDA_CB/CB_gene_counts.csv', row.names = TRUE)

# RNA-seq counts distribution
ggplot(CB_gene_counts) +
  geom_histogram(aes(x = CB0007), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

# Raw counts range
summary(CB_gene_counts)

# few outliers affect distribution visualization
boxplot(CB_gene_counts, main='Raw counts', las=2)

# Rows with few counts, if can be ignored
rs <- rowSums(assay(CB_gse))
hist(log10(rs + 1))
abline(v=1, col="blue", lwd=3)

cts <- CB_gene_counts[rs > 100,]

# two samples’ GE values against each other in a scatterplot
smoothScatter(cts[,1:2]) # better job with so many data points

# take the log of counts & plot again
logcts <- log10(cts + 1)
smoothScatter(logcts[,1:2])

# get a sense for the distribution of counts (log-scale)
hist(logcts[,1])

avg_expr <- rowMeans(CB_gene_counts)

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
logcounts <- log2(CB_gene_counts + 1)

# make a colour vector
statusCols <- case_when(CB_gse$bw=="low" ~ "red3", 
                        CB_gse$bw=="normal" ~ "orange")
# Check distributions of samples using boxplots
boxplot(logcounts,
        xlab="",
        ylab="Log2(Counts)",
        las=2,
        col=statusCols,
        main="Log2(Counts)")
legend("topright", legend = c("Low", "Normal"),fill = c("red3", "orange"), cex=0.7)

logcounts <- log(CB_gene_counts[,3],10) 
boxplot(logcounts, main="", xlab="", ylab="Raw read counts per gene (log10)",axes=FALSE)
d <- density(logcounts)
plot(d,xlim=c(1,8),main="",ylim=c(0,.45),xlab="Raw read counts per gene (log10)", ylab="Density")

# Subset the read counts for the 100 highly expressed genes ----
top.exp <- order(rowMeans(CB_gene_counts), decreasing=TRUE)[1:100]
top.exp <- CB_gene_counts[top.exp,]

# Annotate the samples in the subset with their age (check order with design!)
colnames(top.exp)<- Meta_CB$gestational_age
head(top.exp)
pheatmap(top.exp, cluster_cols = F, labels_row = F)

# Clustering genes most variable across samples ----
VarGenes <- apply(CB_gene_counts, 1, var) # compute variance of each gene

#sort the results by variance in decreasing order and select the top 100 genes
topVarGenes <- names(VarGenes[order(VarGenes, decreasing = T)][1:100])

pheatmap(CB_gene_counts[topVarGenes,], scale = 'row', show_rownames = FALSE)

# overlay some annotation tracks to observer replicate clustering
heatmapColAnnot <- data.frame(colData(CB_gse)[, c("bw", "sex")])
pheatmap(CB_gene_counts[topVarGenes,], scale = 'row',
         show_rownames = FALSE,
         annotation_col = heatmapColAnnot)

## Create the DESeq dataset object
dds <- DESeqDataSet(CB_gse, design = ~ bw)
dds
# variance increases with the average read count
meanSdPlot(assay(dds), ranks = FALSE)

# visualize this design
vd <- VisualizeDesign(sampleData = colData(dds)[, c("sex", "bw")], 
                      designMatrix = attr(dds, "modelMatrix"), 
                      flipCoordFitted = TRUE)
vd

# Transform counts for data visualization ----
CB_rld <- rlog(dds, blind=TRUE)
meanSdPlot(assay(CB_rld), ranks = FALSE)

# Plot PCA, ntop= can be used 
plotPCA(CB_rld, intgroup="bw")

# alternate to show variances
pcaData <- plotPCA(CB_rld, intgroup="sex", returnData=T)

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = sex), size = 3) +
  theme_bw() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() 

pcaData_CB_rld <- DESeq2::plotPCA(CB_rld, intgroup = c("libSize"),
                              returnData = TRUE)
percentVar <- round(100 * attr(pcaData_CB_rld, "percentVar"))

ggplot(pcaData_CB_rld, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = libSize / 1e6), size = 3) +
  theme_bw() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()

# Explore the additional PCs in your data or if you would like to identify genes that contribute most to the PCs. Input is a matrix of log transformed values
CB_rld_mat <- assay(CB_rld)
pcData <- prcomp(t(CB_rld_mat))

# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(Meta_CB, pcData$x)
ggplot(df) +
  geom_point(aes(x=PC1, y=PC3, color = bw)) + theme_bw()

# PCA plot with sample lables 
autoplot(pcData, data = Meta_CB, colour = 'bw', shape='sex', size=3) +
  geom_text_repel(aes(x = PC1, y = PC2, label = cordBlood_ID), box.padding = 0.8)

# count density for all samples, by group
CB_rld_mat %>% 
  as.data.frame() %>% 
  pivot_longer(names_to = "cordBlood_ID", values_to = "logCounts", everything()) %>% 
  left_join(Meta_CB) %>% 
  ggplot(aes(x=logCounts, group = cordBlood_ID)) +
  geom_density(aes(colour = bw)) +
  #scale_colour_manual(values = statusCols) +
  labs(x = "log of Counts", title = "Gene exp count density")

mat <- assay(rld_CB)
mm <- model.matrix(~bw, colData(rld_CB))
mat <- limma::removeBatchEffect(mat, batch=rld_CB$sex, design=mm)
assay(rld_CB) <- mat
plotPCA(rld_CB, intgroup='bw')

# batches in the PCA plot after vst ----
mat <- assay(rld_CB)
mm <- model.matrix(~bw, colData(rld_CB))
mat <- limma::removeBatchEffect(mat, batch=rld_CB$batch, design=mm)
assay(rld_CB) <- mat
plotPCA(rld_CB, intgroup='bw')

# CLUSTERING OF NORMALISED, FILTERED DATA ----
Meta_CB$cl_bw <- as.factor(Meta_CB$bw)
Meta_CB$cl_bw <- factor(Meta_CB$cl_bw, levels=c("normal", "low"))
levels(Meta_CB$cl_bw) <- c("lightblue","darkblue")
Meta_CB$cl_bw <- as.character(Meta_CB$cl_bw)

d <- dist(t(CB_rld_mat))
hc <- hclust(d, method="complete")

# Plot clustering, identifying sample and color coding by status:
pdf("PL/NormalisedFiltered_Clustering_bw.pdf", width=25, height=10)
myplclust(hc, labels=Meta_CB$cordBlood_ID, lab.col=Meta_CB$cl_bw, cex=1.5, main="CordBlood Samples Clustering (complete)")
legend("topright",legend=c("Normal","Low"), cex = 1,
       text.col=c("lightblue","darkblue"), pch=rep(16,2),col=c("lightblue","darkblue"))
dev.off()

# Hierarchical Clustering ----
hclDat <-  t(CB_rld_mat) %>%
  dist(method = "euclidean") %>%
  hclust()
ggdendrogram(hclDat, rotate=TRUE, scale.color=statusCols)

# more customized dendrogram
dendro.dat <-as.dendrogram(hclDat) %>% dendro_data()

dendro.dat$labels <- dendro.dat$labels %>%
  left_join(Meta_CB, by = c(label = "cordBlood_ID"))

ggplot(dendro.dat$segment) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_label(data = dendro.dat$labels,
             aes(x = x,
                 y = y,
                 label = label,
                 fill = bw),
             hjust = 0,
             nudge_y = 1) +
  #scale_fill_manual(values = samgrpCols) +
  coord_flip() +
  labs(x = NULL, y = "Distance", title = NULL, fill = "Sample Group") +
  scale_y_reverse(expand = c(0.3, 0)) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank())
# CLUSTER PLOTS 
bw.bar <- c("grey75","darkblue")
bw.bar <- bw.bar[as.numeric(as.factor(Meta_CB$bw))]

d <- dist(t(CB_rld_mat))
hc <- hclust(d, method="average")
dend <- as.dendrogram(hc)
dend <- color_branches(dend, k=2, col=c("darkgreen","darkred"))

# Draw dendrogram
plot(dend)
colored_bars(colors=bw.bar, dend=dend, rowLabels=c("bw")) # needed ?
legend("topright", legend=c("normal","low"), pch=15, bty="n", col=c("darkblue","grey75"))

# CHECK GENDER CLUSTERING ----

# X: XIST (ENSG00000229807), Y: RPS4Y1 (ENSG00000129824), Y: EIF1AY (ENSG00000198692)
# Y: DDX3Y (ENSG00000067048), Y: KDM5D (ENSG00000012817)

# Subset the normalised counts to the sex-linked genes above:
sex.genes <- c("ENSG00000229807","ENSG00000129824","ENSG00000198692","ENSG00000067048","ENSG00000012817")
norm.sex <- subset(CB_normalized_counts, rownames(CB_normalized_counts) %in% sex.genes)

dim(norm.sex)

# Cluster plot of sex-linked genes:
BW <- c("skyblue","grey75")
BW <- BW[as.numeric(as.factor(Meta_CB$bw))]
Sex <- c("darkred","darkblue")
Sex <- Sex[as.numeric(as.factor(Meta_CB$newborn_gender))]

d <- dist(t(norm.sex))
hc <- hclust(d, method="average")
dend <- as.dendrogram(hc)

pdf("EDA_CB/SexLinkGenes_EuclideanAvgClustering.pdf", width=50, height=10)
par(mar=c(10,7,2,4))
plot(dend)
colored_bars(colors=cbind(BW,Sex), dend=dend, rowLabels=c("BW","Sex"))
dev.off()

# Dendrogram + Heatmap to visualize Euclidean distances ----
raw.dist <- dist(t(assay(CB_rld)))

colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)

ComplexHeatmap::Heatmap(
  as.matrix(raw.dist), 
  col = colors,
  name = "Euclidean\ndistance",
  cluster_rows = hclust(raw.dist),
  cluster_columns = hclust(raw.dist),
  bottom_annotation = columnAnnotation(
    birthweight = CB_rld$bw ,
    sex = CB_rld$sex,
    col = list(sex = c(Female = "pink2", Male = "lightblue3"),
               birthweight  = c(low = "forestgreen", normal = "purple3")))
) 
# rearrange similar samples if clustering not required

# Correlation plot ----
CorMat <- cor(CB_rld_mat) # create cor matrix

corrplot(CorMat, order = 'hclust',
         addrect = 2, addCoef.col = 'white',  # split clusters into group & surround with rectangle
         number.cex = 0.7)
corrplot(CorMat, method = 'square', type = 'lower', diag = FALSE)

# clustering + heatmap of correlation matrix
heat.colors <- brewer.pal(6, "Blues")
pheatmap(CorMat, annotation_col = heatmapColAnnot, color = heat.colors,
         border_color=NA, fontsize = 10, fontsize_row = 10, height=20)

# split the clusters into *two* based on the clustering similarity
pheatmap(CorMat,
         annotation_col = heatmapColAnnot,
         cutree_cols = 2 ) # may use cutree_rows

# Corplot of top var genes in all samples
corVar <- cor(CB_rld_mat[topVarGenes,])
rownames(corVar) <- str_c(Meta_CB$bw, "_", Meta_CB$sex)
colnames(corVar) <- str_c(Meta_CB$bw, "_", Meta_CB$sex)
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
mean_cb <- rowMeans(CB_gene_counts)

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
