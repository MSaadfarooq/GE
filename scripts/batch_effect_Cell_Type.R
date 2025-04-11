# CHECK FOR BATCH EFFECTS ----
library(tidyverse)
library(sva)
library(DESeq2)
library(rafalib)

# load counts and coldata
PL_gene_counts <- read_csv("Outputs/PL_gene_counts.csv")
PL_gene_counts <- column_to_rownames(PL_gene_counts, "geneID")

## Create a matrix of the get normalized counts
dds <- DESeqDataSetFromMatrix(PL_gene_counts, PL_meta, ~bw)
dds$batch <- factor(dds$batch)
table(dds$bw, dds$batch)

dds <- estimateSizeFactors(dds)
norm.cts <- counts(dds, normalized=TRUE)

# Create models (mod=full model, mod0=null, NB, including known batch effects: bw, tissue):
mod <- model.matrix(~as.factor(bw), data=PL_meta)
mod0 <- model.matrix(~1, data=PL_meta)

# Calculate the first 2 SVs (norm.cts - as.matrix())
svseq <- svaseq(norm.cts, mod, mod0, n.sv=2)

# Run DESeq2 correcting for 2 SVs:
dds.2 <- dds
dds.2$SV1 <- svseq$sv[,1]
dds.2$SV2 <- svseq$sv[,2]

# Create Euclidean distance plots for uncorrected & corrected data:
qc.0 <- vst(dds,blind=F)
qc.2 <- vst(dds.2,blind=F)

dist.0 <- dist(t(assay(qc.0)))
dist.2 <- dist(t(assay(qc.2)))

cols <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)

pdf("EDA_PL/SVA_0SVs_EuclideanDistance.pdf", width=25, height=25)
matrix <- as.matrix(dist.0)
rownames(matrix) <- qc.0$bw
colnames(matrix) <- NULL
pheatmap(matrix, clustering_distance_rows=dist.0, clustering_distance_cols=dist.0, col=cols)
dev.off()

pdf("EDA_PL/SVA_2SVs_EuclideanDistance.pdf", width=25, height=25)
matrix <- as.matrix(dist.2)
rownames(matrix) <- qc.2$bw
colnames(matrix) <- NULL
pheatmap(matrix, clustering_distance_rows=dist.2, clustering_distance_cols=dist.2, col=cols)
dev.off()

# Plot the factors estimated by SVA for 2 SVs:

pdf("EDA_PL/SVA_2SVs_BW_Factors.pdf")
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(svseq$sv[,i] ~ dds$bw, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
}
dev.off()

# plot the estimated surrogate variables
bigpar()
dds$bw.int <- as.integer(dds$bw) + 15
plot(svseq$sv[,1:2], col=dds$batch, pch=dds$bw.int, cex=1.5,
     xlab="SV1", ylab="SV2")
legend("top", levels(dds$batch), pch=16,
       col=1:3, cex=.8, ncol=3, title="batch")

# https://biodatascience.github.io/compbio/dist/sva.html
# https://github.com/quadbio/RNAseq_tutorial/

# cell-type deconvulation
library(immunedeconv)
library("AnnotationHub")

ah <- AnnotationHub()
human_ens <- ah[["AH116860"]] #Load the the ensembldb object

# Extract a table of gene annotations
anno.genes <- genes(human_ens,return.type = "data.frame") %>%
  dplyr::select(gene_id,gene_name,gene_biotype) %>% as_tibble()
anno.genes

# merge our result tables with the annotations to add extra columns
CB_counts_anno <- CB_counts %>%
  left_join(anno.genes, by = c("X" = "gene_id"))

CB_counts_hgnc <- CB_counts_anno |> dplyr::select(-c("ensid","X","gene_biotype"))

# Handle Duplicate Gene Names
CB_counts_hgnc <- CB_counts_hgnc %>%
  distinct(gene_name, .keep_all = TRUE) %>% # Keep only first occurrence of each gene
  column_to_rownames(var = "gene_name")

deconvolution_methods

# obtain a cell_type x sample data frame with cell-type scores for each sample
res_quantiseq <- deconvolute(gene_expression = CB_counts_hgnc, method = "quantiseq")

# QuanTIseq generates scores that can be interpreted as a cell-type fraction. Let’s visualize
res_quantiseq |>
  pivot_longer(cols = -cell_type, names_to = "sample", values_to = "fraction") |> 
  ggplot(aes(x = sample, y = fraction, fill = cell_type)) +
  geom_bar(stat = "identity") + # stacked bar-chart suggests the scores to be cell-type fractions
  coord_flip() +
  scale_fill_brewer(palette = "Paired") +
  scale_x_discrete(limits = rev(levels(res_quantiseq)))

# MCP-counter: scores in arbitrary units, only comparable between samples, but not between cell-types
res_mcp_counter <- deconvolute(CB_counts_hgnc, "mcp_counter")

# visualize the scores per-cell type, allowing for a relative comparison between samples
res_mcp_counter %>%
  pivot_longer(cols = -cell_type, names_to = "sample", values_to = "score") %>%
  ggplot(aes(x = sample, y = score, color = cell_type)) +
  geom_point(size = 2) +
  facet_wrap(~cell_type, scales = "free_x", ncol = 3) +
  scale_color_brewer(palette = "Paired", guide = FALSE) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

set_cibersort_binary("/path/to/CIBERSORT.R")
set_cibersort_mat("/path/to/LM22.txt")

deconvolute(CB_counts_hgnc, "cibersort")   # or 'cibersort_abs'