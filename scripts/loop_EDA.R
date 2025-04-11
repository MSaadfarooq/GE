# EDA Loop

# Load necessary libraries
library(pheatmap)

# Read the metadata file
meta <- read.csv("~/rna_seq/Meta_CB.csv")

# Get a list of all gene count files (assuming they are named like 'gene_counts_*.csv')
file_list <- list.files(path = "~/rna_seq", pattern = "*\\_gene_counts.csv$", full.names = TRUE)

# Loop over each gene count file
for (file in file_list) {
  # Extract the tissue name from the file name (without path and extension)
  tissue_name <- gsub("*\\_gene_counts.csv$", "\\1", basename(file))
  
  # Read the current gene count file
  gene_counts <- read.csv(file, row.names = 1)
  
  # Print summary of gene counts (for debugging or checking)
  summary(gene_counts)
  
  # Create a boxplot with the tissue name in the title
  boxplot(gene_counts, main = paste('Raw counts for', tissue_name), las = 2)
  
  # Extract the top 100 most expressed genes
  top.exp <- order(rowMeans(gene_counts), decreasing = TRUE)[1:100]
  top.exp <- gene_counts[top.exp,]
  
  # Saving the heatmap to a PDF
  pdf(paste("heatmap_", tissue_name, ".pdf", sep = ""))
  pheatmap(top.exp, cluster_cols = FALSE, labels_row = FALSE, main = paste('Heatmap for', tissue_name))
  dev.off()
}
###################
# Load necessary libraries
library(DESeq2)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)

# Get a list of all 'gse' object RDS files (assuming they are saved as 'gse_*.RDS')
gse_files <- list.files(path = "~/rna_seq/Outputs", pattern = "gse_.*\\.RDS$", full.names = TRUE)

# Loop over each 'gse' object RDS file
for (gse_file in gse_files) {
  
  # Extract the tissue name from the file name (without path and extension)
  tissue_name <- gsub("gse_(.*)\\.RDS", "\\1", basename(gse_file))
  
  # Load the SummarizedExperiment (gse) object
  gse <- readRDS(gse_file)
  
  # Compare the library sizes of samples
  gse$libSize <- colSums(assay(gse))
  
  # Create a library size plot
  ggplot(as.data.frame(colData(gse)), aes(x = rownames(colData(gse)), y = libSize / 1e6, fill = bw)) + 
    geom_bar(stat = "identity") + theme_bw() + 
    labs(x = "Sample", y = "Total count in millions") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  # Write the counts to an object
  gene_counts <- assay(gse) %>% 
    round() %>% 
    data.frame()
  
  # Save the current counts and gse object
  write.csv(gene_counts, paste0('Outputs/gene_counts_', tissue_name, '.csv'), row.names = TRUE)
  
  # Create a boxplot with the tissue name in the title
  boxplot(gene_counts, main = paste('Raw counts for', tissue_name), las = 2)
  
  # Extract the top 100 most expressed genes
  top.exp <- order(rowMeans(gene_counts), decreasing = TRUE)[1:100]
  top.exp <- gene_counts[top.exp,]
  
  # Saving the heatmap to a PDF
  pdf(paste("heatmap_", tissue_name, ".pdf", sep = ""))
  pheatmap(top.exp, cluster_cols = FALSE, labels_row = FALSE, main = paste('Heatmap for', tissue_name))
  dev.off()
  
  # Log2 transformation of counts
  logcounts <- log2(gene_counts + 1)
  
  # Make a color vector for boxplot
  statusCols <- case_when(gse$bw == "low" ~ "red3", 
                          gse$bw == "normal" ~ "orange")
  
  # Boxplot for log-transformed counts
  boxplot(logcounts,
          xlab = "",
          ylab = "Log2(Counts)",
          las = 2,
          col = statusCols,
          main = paste("Log2(Counts) for", tissue_name))
  
  legend("topright", legend = c("Low", "Normal"), fill = c("red3", "orange"), cex = 0.7)
  
  # Create DESeqDataSet object
  dds <- DESeqDataSet(gse, design = ~ bw)
  
  # Transform counts for data visualization using VST
  vsd <- vst(dds, blind = TRUE)
  
  # Save transformed data (VST)
  saveRDS(vsd, file = paste0("Outputs/vsd_", tissue_name, ".RDS"))
  
  # Mean-Standard Deviation plot
  meanSdPlot(assay(vsd), ranks = FALSE)
  
  # PCA plot
  pcaData <- plotPCA(vsd, intgroup = "bw", returnData = TRUE)
  
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  # PCA plot with percent variance
  PC_1_2 <- ggplot(pcaData, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = bw), size = 3) +
    theme_bw() +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed()
  ggsave(filename = paste("figures/", tissue_name, "_PC_var.png", sep = ""), plot = PC_1_2)
  
  vsd_mat <- assay(vsd)
  
  pcData_add <- prcomp(t(vsd_mat))
  
  # Create data frame with metadata and additional PC values for input to ggplot
  df <- cbind(meta, pcData_add$x)
  PC_add <- ggplot(df) +
    geom_point(aes(x=PC1, y=PC3, color = group)) + theme_bw()
  ggsave(filename = paste("figures/", tissue_name, "_PC_add.png", sep = ""), plot = PC_add)
  
  # CLUSTERING OF NORMALISED, FILTERED DATA ----
  meta$cl <- as.factor(meta$bw)
  meta$cl <- factor(meta$cl, levels=c("normal", "low"))
  levels(meta$cl) <- c("lightblue","darkblue")
  meta$cl <- as.character(meta$cl)
  
  d <- dist(t(vsd_mat))
  hc <- hclust(d, method="complete")
  
  # Plot clustering, identifying sample and color coding by status
  pdf(paste("NormalisedFiltered_Clustering_", tissue_name, ".pdf", sep = "")) #width=25, height=10
  myplclust(hc, labels=tissue_name, lab.col=meta$cl, cex=1.5, main="Samples Clustering (complete)")
  legend("topright",legend=c("normal", "low"), cex = 1,
         text.col=c("lightblue","darkblue"), pch=rep(16,2),col=c("lightblue","darkblue"))
  dev.off()
  
  hclDat <-  t(vsd_mat) %>% dist(method = "euclidean") %>%
    hclust()
  pdf(paste("Norm_Filt_Clustering_", tissue_name, ".pdf", sep = ""))
  ggdendrogram(hclDat, rotate=TRUE, scale.color=statusCols)
  dev.off()
  
  # more customized dendrogram
  dendro.dat <- as.dendrogram(hclDat) %>% dendro_data()
  
  dendro.dat$labels <- dendro.dat$labels %>%
    left_join(meta, by = c(label = "sample"))
  
  dendro.custom <- ggplot(dendro.dat$segment) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_label(data = dendro.dat$labels,
               aes(x = x, y = y, label = label,fill = group),hjust = 0, nudge_y = 1) +
    #scale_fill_manual(values = samgrpCols) +
    coord_flip() +
    labs(x = NULL, y = "Distance", title = NULL, fill = "Sample Group") +
    scale_y_reverse(expand = c(0.3, 0)) +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), panel.background = element_blank())
  ggsave(filename = paste("figures/", tissue_name, "_dendro.custom.png", sep = ""), plot = dendro.custom)
  
  # Dendrogram + Heatmap to visualize Euclidean distances ----
  raw.dist <- dist(t(assay(vsd)))
  colz <- colorRampPalette(brewer.pal(9, "Blues"))(255)
  pdf(paste("Heatmap_Cluster_", tissue_name, ".pdf", sep = ""))
  Heatmap(as.matrix(raw.dist), col = colz,
          name = "Euclidean\ndistance",
          cluster_rows = hclust(raw.dist),
          cluster_columns = hclust(raw.dist),
          bottom_annotation = columnAnnotation(birthweight = vsd$bw,
                                               col = list(sex = c(Female = "pink2", Male = "lightblue3"),
                                                          birthweight  = c(low = "forestgreen", normal = "lightgreen")))
  )
  dev.off()
  CorMat <- cor(vsd_mat) # create cor matrix
  pdf(paste("corrPlot_", tissue_name, ".pdf", sep = ""))
  corrplot(CorMat, order = 'hclust',
           addrect = 2, addCoef.col = 'white',  # split clusters into group & surround with rectangle
           number.cex = 0.7)
  dev.off()
}
