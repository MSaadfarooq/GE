# Load libraries ----
suppressPackageStartupMessages({
  library(DESeq2)
  library(vsn)
  library(tidyverse)
  library(apeglm)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(ggvenn)
})
## Simple project ID
project="EDA_PL"

# Re-create DESeq2 dataset if the design formula has changed after QC
dds_PL <- DESeq2::DESeqDataSet(PL_gse, design = ~ sex + bw)
dds_PL
## Run DESeq ----
dds_PL <- DESeq(dds_PL)

# Plot dispersions
plotDispEsts(dds_PL, main="Dispersion plot")

# approach-to-count-outliers
pdf("RAW_COUNTS_QC/GSE57945_CooksCutoff_NoBatchCorrection.pdf", width=40, height=7)
boxplot(log10(assays(dds_PL)[["cooks"]]), range=0, las=2)
dev.off()

# dds object transformation
rld_PL <- rlog(dds_PL, blind=F)
saveRDS(rld_PL,file =paste0(project,".rld.RDS"))
saveRDS(dds_PL,file =paste0(project,".dds.RDS"))

# Extract DESeq2 results
resultsNames(dds_PL)
PL_result <- results(dds_PL, alpha=0.05)

summary(PL_result, alpha=0.05)

head(PL_result[order(PL_result$pvalue), ]) # order by lowest p-value

# Output normalized counts to save as a file
PL_normalized_counts <- counts(dds_PL, normalized=TRUE)
write.csv(PL_normalized_counts, 'Placenta_normalized_counts.csv')

# filtering threshold that has been used to filter low count genes
metadata(PL_result)$filterThreshold #genes with basemean < [output] have been filtered

as_tibble(metadata(PL_result)$filterNumRej) %>%
  ggplot(aes(x = theta, y = numRej)) +
  geom_point() +
  geom_vline(xintercept = 0.232,
             color = 'red')

# distribution of adjusted p-values
hist(PL_result$padj, col="lightblue", main = "Adjusted p-value distribution")

# distribution of p-values
hist(PL_result$pvalue, col="royalblue", main = "Non-adjusted p-value distribution")

## Gene-level filtering ----
# Filter genes by zero expression
PL_result[which(PL_result$baseMean == 0),] %>% 
  data.frame() %>% 
  View()

# Filter genes that have an extreme outlier - needs replicate
PL_result[which(is.na(PL_result$pvalue) & 
                  is.na(PL_result$padj) &
                  PL_result$baseMean > 0),] %>% 
  data.frame() %>% View()

# explore the relationship between LFC, significant DEGs and the genes mean count
plotMA(PL_result)

# Independent Filtering and log-fold shrinkage
PL_result_Lfc <- lfcShrink(dds_PL, coef = "bw_normal_vs_low", type="apeglm")
plotMA(PL_result_Lfc)

head(as.data.frame(PL_result))
head(as.data.frame(rowRanges(PL_gse)))

## Extract DE genes ----
# Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 1
PL_DElist <- PL_result %>%
  data.frame() %>%
  tibble::rownames_to_column(var="gene") %>% 
  mutate(threshold=case_when(padj<=padj.cutoff & abs(log2FoldChange)>=lfc.cutoff~"Significant",T~"Not Significant")) %>%
  arrange(padj) %>% as_tibble()

# Subset the tibble to keep only significant genes
sigDE_PL <- PL_DElist %>%
  dplyr::filter(threshold=="Significant")

# or to get df of p-adj
sig_DE_df <- as.data.frame(subset(PL_result, padj < 0.05))
sig_DE_genes <- rownames(sig_DE_df)
length(sig_DE_genes)

# add annotations like the gene name and biotype ----
library(AnnotationHub)
ah <- AnnotationHub() # create a link to the annotation hub
query(ah,c("EnsDb","112","Homo")) # Search annotationHub for Ensembl annotation

human_ens <- ah[["AH116860"]] #Load the the ensembldb object

# Extract a table of gene annotations 
anno.genes <- genes(human_ens,return.type = "data.frame") %>% 
  dplyr::select(gene_id,gene_name,gene_biotype,entrezid) %>% 
  as_tibble()
anno.genes

# resolve entrez-Ensembl ID conflict by selecting the first in each list
anno.genes$entrezid <- map(anno.genes$entrezid,1) %>% unlist()
anno.genes

# merge our result tables [sig or all DE] with the annotations to add extra columns
sigDE_anno <- sig_DE_res %>%
  left_join(anno.genes, by = c("gene" = "gene_id")) %>%
  mutate(gene_biotype = as.factor(gene_biotype)) # biotype as a factor

write.csv(sigDE_anno, file = "PL.sig.DEGs.csv")
readr::rds(sigDE_anno, 'PL_sigDE_anno.RDS')

##===========================##
# Visualization of DE results
## ==========================##

## Volcano plot ----
PL_DElist %>%
  dplyr::filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = threshold)) +
  scale_colour_manual(values = c("gray", "red")) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) + # lines optional
  theme_bw()

## Barplot of gene_biotype
ggplot(sigDE_anno,aes("Significant genes",fill=gene_biotype)) +
  geom_bar() + 
  theme_bw() + 
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## Heatmap of results ----
# Transform counts
rld_PL <- rlog(dds, blind=F)

# Get top 300 DE genes
HM_genes <- PL_result[order(PL_result$padj), ] |>
  head(300) |>
  rownames()
heatmapData <- assay(rld_PL)[HM_genes, ]

# Scale counts for visualization
heatmapData <- t(scale(t(heatmapData)))

# Add annotation & colors
heatmapColAnnot <- data.frame(colData(rld_PL)[, c("bw", "sex")])
heatmapColAnnot <- HeatmapAnnotation(df = heatmapColAnnot)
colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)

# Plot as Heatmap
ComplexHeatmap::Heatmap(heatmapData,
                        col = colors,
                        top_annotation = heatmapColAnnot,
                        cluster_rows = TRUE, cluster_columns = FALSE)

# Top 10 differentially expressed genes (by padj values) ----
# Order results by padj values
top10_sigOE <- PL_DElist %>% 
  arrange(padj) %>% 	# Arrange rows by padj
  #pull(gene) %>% 		# Extract character vector of ordered genes
  head(10)		# Extract the first 10 genes

meta<-coldata %>% rename(sample=names)
  
as_tibble(counts(dds[top10_sigOE$gene, ]),
          rownames = 'gene') %>%
  pivot_longer(names_to = "sample", values_to = "counts", -gene) %>%  
  left_join(as_tibble(coldata, rownames = "sample")) %>% 
  ggplot(aes(x = sample, y = counts, fill = bw)) + #geom_point()
  geom_bar(stat = 'identity', color = "gray30") +
  facet_wrap( ~ gene, scales = "free", ncol = 3) +
  theme(axis.text.x = element_text(size = 7, angle = 90),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7))

